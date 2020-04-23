#include <algorithm>
#include <memory>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <iomanip>
#include <complex>
#include <cassert>
#include <utility>
#include "advection/example_meshes.h"
#include "advection/transport.h"
#include <mpi.h>
#include <omp.h>
#include <cuda_runtime.h>
#include <ctime>
#include <cstdlib>

#include "helper.cuh"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Triangulation_2.h>
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
typedef CGAL::Triangulation_2<K,Tds> Triangulation;
typedef CGAL::Triangle_2<K> Triangle;
typedef CDT::Face Face;
typedef CDT::Face_handle Face_handle;
typedef CDT::Finite_faces_iterator Face_iterator;


using namespace std;
using namespace example_meshes;
using namespace advection;

int get_distance(CDT& cdt, Face_handle face)
{
	int count = 0;
	for (auto itt(cdt.finite_faces_begin()); itt != cdt.finite_faces_end(); itt++)
	{
		Face_handle dd = itt;
		if (dd==face)
		{
			return count;
		}
		count++;
	}
	std::cout << "IT'S NOT FOUND" << std::endl;
	return count;
}

Face_handle get_find(CDT& cdt, Face_handle face)
{
	int count = 0;
	for (auto itt(cdt.finite_faces_begin()); itt != cdt.finite_faces_end(); itt++)
	{
		Face_handle dd = itt;
		if (dd==face)
		{
			return itt;
		}
		count++;
	}
	std::cout << "IT'S NOT FOUND" << std::endl;
	return nullptr;
}

int main(int argc, char ** argv)
{
	MPI_Init(&argc, &argv);
	int size_of_group;
    	int rank;
	MPI_Comm comm = MPI_COMM_WORLD;
    	MPI_Comm_size(comm, &size_of_group);
    	MPI_Comm_rank(comm, &rank);
	if (rank==0) cout << "mpi 2d version" << endl;

	double h, dr, dt;
	int TIME, TOTAL_FRAMES;
	double length_x, length_y;
	double epsilon_1, epsilon_2, mu_1, mu_2;

	int N_CUDA;
	int N_THREADS_PER_PROC;
	
	if(argc !=2) {cout << "need filename" << endl; return 2;}

	string filename{argv[1]};
	ifstream arguments;
	arguments.open(filename);
	vector<string> data;
	string line;
	int i = 0;
	while(getline(arguments, line))
	{
		data.push_back(line);
		i++;
	}
	if (data.size() != 13) {cout << "wrong data" << endl; return 2;}
	h = stod(data[0]);//0.05
	dr = stod(data[1]);//0.04
	dt = stod(data[2]);//0.001
	TIME = stoi(data[3]);//1
	TOTAL_FRAMES = stoi(data[4]);//100
	length_x = stod(data[5]);//2.0
	length_y = stod(data[6]);//0.5
	epsilon_1 = stod(data[7]);
	epsilon_2 = stod(data[8]);
	mu_1 = stod(data[9]);
	mu_2 = stod(data[10]);
	N_CUDA = stoi(data[11]);
	N_THREADS_PER_PROC = stoi(data[12]);

	bool if_paint = (TOTAL_FRAMES == 0) ? 0 : 1;
	int periods = (TOTAL_FRAMES == 0) ? TIME : TIME/TOTAL_FRAMES;

	
	if (rank >= N_CUDA || N_CUDA==0)
	{
		omp_set_num_threads(N_THREADS_PER_PROC);
	}
	if (rank < N_CUDA) set_device(rank);
	double * proportions = new double[4*size_of_group];
	
	CDT cdt;
	vector<Face_handle> visual_core;
	for (int i = 0; i < size_of_group; i++) 
	{
		double x_length, y_length;
		if (N_CUDA==0 || size_of_group - N_CUDA == 0)
		{
			x_length = 1.5;
			y_length = 0.7;
		}
		else if (i < N_CUDA)
		{
			x_length = 1.2;
			y_length = 0.7;
		}
		else
		{
			x_length = 0.3;
			y_length = 0.7;
		}
		double x_len = x_length;
		double y_len = y_length;
		int x_times = 0;
		int y_times = 0;
		if (i<N_CUDA)
		{
			while ( (y_length/y_len)*(x_length/x_len) != N_CUDA)
			{
				if (x_len > y_len)
				{
					x_len/=2.0;
					x_times++;
				}
				else
				{
					y_len/=2.0;
					y_times++;
				}
			}

			proportions[4*i] = x_times!=0 ? x_len*(i%(x_times*2)) : 0.0;
			proportions[4*i+1] = x_times!=0 ? x_len+x_len*(i%(x_times*2)) : x_len;
			proportions[4*i+2] =x_times!=0 ?( y_len* ((int) (i/(x_times*2))) ): (y_len * i);  //y_len*(i/(size_of_group/2/y_times/2) );
			proportions[4*i+3] = x_times!=0 ? ( y_len* ((int) (i/(x_times*2))) + y_len): (y_len + y_len * i);
		}
		else
		{
			while ( (y_length/y_len)*(x_length/x_len) != (size_of_group - N_CUDA) )
			{
				if (x_len > y_len)
				{
					x_len/=2.0;
					x_times++;
				}
				else
				{
					y_len/=2.0;
					y_times++;
				}
			}
			proportions[4*i] = (x_times!=0 ? x_len*(i%(x_times*2)) : 0.0) + 1.5 - x_length;
			proportions[4*i+1] = (x_times!=0 ? x_len+x_len*(i%(x_times*2)) : x_len) + 1.5 - x_length;
			proportions[4*i+2] = x_times!=0 ? ( y_len* ((int) ((i-N_CUDA)/(x_times*2))) ): (y_len * (i-N_CUDA) );  //y_len*(i/(size_of_group/2/y_times/2) );
			proportions[4*i+3] = x_times!=0 ? ( y_len* ((int) ((i-N_CUDA)/(x_times*2))) +y_len ): (y_len * (i-N_CUDA) + y_len);
		}
	}
	if (rank == 0) for (int i = 0; i < size_of_group; i++) cout << "here " << proportions[4*i+0] << " " << proportions[4*i+1] << " " << proportions[4*i+2] << " " << proportions[4*i+3] << endl;

	int ** ids_to_send = new int* [size_of_group];
	int * ids_size = new int [size_of_group];
	mesh_mpi(cdt, dr, length_x, length_y, size_of_group, proportions, visual_core, ids_to_send, ids_size, rank);
	//cout << "rank is " << rank << " number of triangles is " << visual_core.size() << endl;

	if (rank==0 && if_paint)
	{
		std::ofstream verts;
  		verts.open("vertices.txt");
  		std::ofstream faces;
  		faces.open("faces.txt");
	
  		for (auto it(cdt.points_begin()); it != cdt.points_end(); it++)
  		{
			verts << (it)->x() << " " << (it)->y() << " 0.0" << std::endl; 
  		}
  		verts.close();
  
  		for (auto it(cdt.finite_faces_begin()); it != cdt.finite_faces_end(); it++)
  		{
			for (int i = 0; i < 3; i++)
			{
				auto vert = std::find(cdt.points_begin(), cdt.points_end(), cdt.triangle(it)[i]);
				if (vert != cdt.points_end()) faces << std::distance(cdt.points_begin(), vert) << " ";
				else std::cout << "SOMETHING HAPPENED" << std::endl;
			}
			faces << std::endl;
  		}
  		faces.close();
	}
	//cout << "finished drawing vertices and face ids" << endl;
	int ** ids;
	if (rank==0)
	{ 
		int *buffer;
		ids = new int * [size_of_group];
        	for (int i = size_of_group-1; i >= 0; i--)
		{
			if (i != 0 )
			{
				MPI_Status status;
				MPI_Probe(i, 1, comm, &status);
				int size = 0;
				MPI_Get_count(&status, MPI_INT, &size);
				buffer = new int [size];

				ids[i] = new int [size];
				MPI_Recv(buffer, size, MPI_INT, i, 1, comm, &status);
				for (int k = 0; k < size; k++)
				{	
					ids[i][k] = buffer[k];
				}
				delete[] buffer;
			}
			else
			{
				ids[i] = new int [visual_core.size()];
				int k = 0;
				for (auto it(visual_core.begin()); it !=visual_core.end(); it++)
				{
					auto global_tr = get_find(cdt, *it);
					ids[i][k] = get_distance(cdt, global_tr);
					k++;
				}
			}
		}
    	} 
	else 
	{
		int * id = new int[visual_core.size()];
		int k = 0;
		for (auto it(visual_core.begin()); it !=visual_core.end(); it++)
		{
			auto global_tr = get_find(cdt, *it);
			id[k] = get_distance(cdt, global_tr);
			k++;
		}
        	MPI_Ssend(id, visual_core.size(), MPI_INT, 0, 1, comm);
    	}

	//DATA FOR PYTHON VISUALIZER BEGIN
	ofstream colours;
	if (rank==0 && if_paint) colours.open("colours.txt", std::ios_base::app);
	advection_2d_mpi * equation;
	equation = new advection_2d_mpi(rank<N_CUDA ? 1:0, rank, size_of_group, cdt, visual_core, ids_to_send, ids_size);	

	double *E_y;	
	double *E_x;
	double *H_z;
	double * coefs_E_y;
	double * coefs_E_x;
	double * coefs_H_z;
	double *epsilon;
	double *mu;
	
	if (rank < N_CUDA)
	{
		initialize_cuda_memory(E_y, E_x, H_z, coefs_E_x, coefs_E_y, coefs_H_z, mu, epsilon, visual_core.size());
	}
	else
	{
		E_y = new double [visual_core.size()];
		E_x = new double [visual_core.size()];
		H_z = new double [visual_core.size()];
		coefs_E_y = new double [2*visual_core.size()];
		coefs_E_x = new double [2*visual_core.size()];
		coefs_H_z = new double [2*visual_core.size()];
		epsilon = new double [visual_core.size()];
		mu = new double [visual_core.size()];
	}	
	
	double *u_collect;
	if (rank == 0)
	{
		u_collect = new double [cdt.number_of_faces()];
	}

	int p = 0;
	for (auto it(visual_core.begin()); it != visual_core.end(); it++)
	{
		E_x[p] = 0.0;
		E_y[p] = 0.0;
		H_z[p] = 0.0;
		epsilon[p] = pow(CGAL::centroid(cdt.triangle(*it)).x()- 0.75,2)+pow(CGAL::centroid(cdt.triangle(*it)).y()- 0.35,2) < 0.01 ? epsilon_2 : epsilon_1;
		mu[p] = pow(CGAL::centroid(cdt.triangle(*it)).x()- 0.75,2)+pow(CGAL::centroid(cdt.triangle(*it)).y()- 0.35,2) < 0.01 ? mu_2 : mu_1;
		p++;
	}	
	
	double start;
	MPI_Barrier(comm);
	if(rank==0) cout << "starting this stuff..." << endl;
	if(rank == 0) start = MPI_Wtime();
	for (int t = 0; t < TIME; t++)
	{
		if(rank < N_CUDA)
		{
			//cout << t << endl;
			call_update_coefs_cuda(H_z, E_x, E_y, coefs_H_z, coefs_E_x, coefs_E_y, mu, epsilon, visual_core.size());
		}
		else
		{
			#pragma omp parallel for
			for (int p = 0; p < visual_core.size(); p++)
			{
				//if (t==0) {int id = omp_get_thread_num();

        			//printf("rank: %d\tthread: %i\n",rank, id);}
				coefs_H_z[2*p] = 1.0/mu[p]*(E_y[p]);
				coefs_H_z[2*p+1] = -1.0/mu[p]*E_x[p];
				
				coefs_E_x[2*p] = 0.0;
				coefs_E_x[2*p+1] = -1.0/epsilon[p]*(H_z[p]);

				coefs_E_y[2*p] = 1.0/epsilon[p]*(H_z[p]);
				coefs_E_y[2*p+1] = 0.0;
			}
		}
		MPI_Barrier(comm);
		equation->solver_small(H_z, dt, coefs_H_z, true, true, t*dt);
		equation->solver_small(E_x, dt, coefs_E_x, true, false, t*dt);
		equation->solver_small(E_y, dt, coefs_E_y, false, false, t*dt);
		if (t % periods == 0 && if_paint)
		{
    			if (rank==0)
			{ 
				double *buffer;
        			for (int i = size_of_group-1; i >= 0; i--)
				{
					if (i != 0 )
					{
						MPI_Status status;
						MPI_Probe(i, 1, comm, &status);
						int size = 0;
						MPI_Get_count(&status, MPI_DOUBLE, &size);
						buffer = new double [size];
						MPI_Recv(buffer, size, MPI_DOUBLE, i, 1, comm, &status);
							//cout << "size " << size << endl;
							for (int k = 0; k < size; k++)
							{	
								u_collect[ids[i][k]] = buffer[k];
							}
						delete[] buffer;
					}
					else
					{
						for(int k = 0; k < visual_core.size(); k++)
							u_collect[ids[i][k]] = H_z[k];
					}
				}
    			} 
			else 
			{
            			MPI_Ssend(H_z, (visual_core.size()), MPI_DOUBLE, 0, 1, comm);
    			}
		}
		if (t%periods == 0 && if_paint && rank==0) 
		{
			double * p_n = &u_collect[0];
			double maximum = 0.0;
			double minimum = 0.0;
			double average = 0.0;
			for (int k = 0; k < cdt.number_of_faces(); k++)
			{
				if (p_n[k] > maximum) maximum = p_n[k];
				if (p_n[k] < minimum) minimum = p_n[k];
				average+=p_n[k];
			}
			average = average/cdt.number_of_faces();
			//cout << minimum << endl;
			//cout << maximum << endl;
			//cout << average << endl;
			double * p_new = new double [cdt.number_of_faces()];
			for (int k = 0; k < cdt.number_of_faces(); k++)
			{
				//p_new[k] = p_n[k];
				p_new[k] = (p_n[k] - -1)/(1 - -1);
			}
			for (int k = 0; k < cdt.number_of_faces(); k++)
			{
				for (int i = 0; i < 1; i++)
				{
					int index = k+cdt.number_of_faces()*i;
					if ((p_new[index]) <= 0.25) colours << 0.0 << " " << ((p_new[index])/0.25)  << " " << 1.0 << " " << 1.0 << " "; 
					else if ((p_new[index]) <= 0.5) colours << 0.0 << " " << 1.0 << " " << ((0.5-(p_new[index]))/0.25) << " " << 1.0 << " ";
					else if ((p_new[index]) <= 0.75) colours << (((p_new[index])-0.5)/0.25) << " " << 1.0 << " " << 0.0 << " " << 1.0 << " ";
					else colours << 1.0 << " " << ((1.0-p_new[index])/0.25) << " " << 0.0 << " " << 1.0 << " ";	
				}
				colours << endl;
			}
			delete [] p_new;
		}
	}
	MPI_Barrier(comm);	
	double end;
	if (rank == 0)
	{
		end = MPI_Wtime();
		cout << "duration " << end-start << endl;
		delete [] u_collect;
	}
	if (rank == 0 && if_paint) colours.close();
	MPI_Finalize();
	return 0;
}

