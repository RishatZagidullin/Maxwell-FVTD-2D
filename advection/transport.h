#pragma once
#include <mpi.h>
#include <math.h>	
#include <cuda_runtime.h>
#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <complex>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <iostream>
#include <fstream>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/centroid.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;
typedef CGAL::Vector_2<K> Vector;
typedef CGAL::Triangle_2<K> Triangle;
typedef CDT::Face Face;

typedef CDT::Face_handle Face_handle;
typedef CDT::Finite_faces_iterator Face_iterator;

using namespace std;

namespace advection
{
	class advection_2d
	{
	protected:
		int solver_type;
		int size;
		double * face_normals;
		double * tr_areas;
		double * xbarys;
		double * ybarys;
		int * neighbor_ids;
		double limiter(double const& r_factor, double const& weight);
		void solver_small_cpu(double *&u, double const &dt, double *&velocities, double *&interpolated_velocities, bool if_y, bool if_h, double t);

	public:
		bool * is_boundary;
		void calc_dynamic_arrays(CDT& cdt, vector<Face_handle>& triangles);
		advection_2d(CDT& cdt, vector<Face_handle>& triangles, int solver_type) : size(triangles.size()), solver_type(solver_type)
		{
			calc_dynamic_arrays(cdt, triangles);
		}
		~advection_2d();
		virtual void solver_small(double *&u, double const dt, double *&coefs, bool if_y, bool if_h, double t);
	};

	class advection_2d_mpi : public advection_2d
	{
	protected:
		int rank;
		int size_of_group;
		int * interfaces_size;

		//these two are determined in calc_mpi_arrays
		bool ** send_to;
		int ** received_ids;

		int ** ids_to_send;

		void calc_mpi_arrays(CDT& cdt, vector<Face_handle>& local_triangles);
	public:
		int * global_ids_reference;
		advection_2d_mpi(int solver_type, int rank, int size_of_group, CDT& cdt, vector<Face_handle> local_triangles, int ** & ids_to_send, int * &ids_size) : advection_2d(cdt, local_triangles, solver_type), rank(rank), interfaces_size(ids_size), size_of_group(size_of_group), ids_to_send(ids_to_send)
		{
			calc_mpi_arrays(cdt, local_triangles);
		}
		~advection_2d_mpi();//DON"T FORGET TO WRITE IT OUT, AND REMEMBER IT DOESN"T EXIST ANYWHERE ELSE
		virtual void solver_small(double *& u, double const& dt, double*& coefs, bool if_y, bool if_h, double t);
	};

}
