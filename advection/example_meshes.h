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
typedef CGAL::Triangle_2<K> Triangle;
typedef CDT::Face Face;
typedef CDT::Face_handle Face_handle;
typedef CDT::Finite_faces_iterator Face_iterator;

using namespace std;

namespace example_meshes
{
	void mark_triangle(CDT& cdt, int * proc_labels, int i, int& val, queue<Face_handle> & faces);
	void mesh_mpi(CDT& GRID, double dr, double x_length, double y_length, int n_processes, double const * proportions, vector<Face_handle>& output, int ** & ids_to_send, int * & ids_size, int rank);
	void write_mesh(CDT& GRID, double dr, double x_length, double y_length, int n_processes, double const * proportions, int ** ids);
}
