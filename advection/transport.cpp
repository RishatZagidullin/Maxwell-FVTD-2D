#include "transport.h"
#include "transport.cuh"
namespace advection
{
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
		//std::cout << "IT'S NOT FOUND" << std::endl;
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
		//std::cout << "IT'S NOT FOUND" << std::endl;
		return nullptr;
	}

	double sigma(double const &r, double const &dr, int const &PML, int const & size)
	{
		if (r>(dr*(size-PML))) return pow((r-(size-PML)*dr)/(PML*dr),2)*3.0*log(10.0)*13.0/(PML*dr);
		else if (r < dr*(PML)) return pow((PML*dr-r)/(PML*dr),2)*3.0*log(10.0)*13.0/(PML*dr);
		else return 0;
	}

	advection_2d::~advection_2d()
	{
		dealloc_arrays(solver_type, is_boundary, face_normals, neighbor_ids, tr_areas, xbarys, ybarys);
	}

	advection_2d_mpi::~advection_2d_mpi()
	{
		for (int i = 0; i < size_of_group; i++)
		{
			delete [] send_to[i];
			delete [] ids_to_send[i];
			delete [] received_ids[i];
		}
		delete [] send_to;
		delete [] received_ids;
		delete [] ids_to_send;
		delete [] interfaces_size;
		delete [] global_ids_reference;
	}

	double advection_2d::limiter(double const &r_factor, double const &weight)
	{
		return max(0.0, min(weight*r_factor,min(0.5*(1.0+r_factor), weight)));//MC
	}
	

	void advection_2d::calc_dynamic_arrays(CDT& cdt, vector<Face_handle>& triangles)
	{
		ofstream barys;
		barys.open("cellcenters.txt");
		initialize_arrays(solver_type, size, is_boundary, face_normals, neighbor_ids, tr_areas, xbarys, ybarys);
		int k = 0;
		double leftmost_x = cdt.triangle(triangles[0])[0].x();
		double rightmost_x = leftmost_x;
		double bottom_y = cdt.triangle(triangles[0])[0].y();
		double top_y = bottom_y;
		int counter = 0;
		for (auto it(triangles.begin()); it != triangles.end(); it++)
		{
			barys << -CGAL::centroid(cdt.triangle(*it)).y() << " " << CGAL::centroid(cdt.triangle(*it)).x() << " " << 0.0 << endl;

			xbarys[counter] = CGAL::centroid(cdt.triangle(*it)).x(); 
			ybarys[counter] = CGAL::centroid(cdt.triangle(*it)).y();
			for (int i = 0; i < 3; i++)
			{
				double value = cdt.triangle(*it)[i].x();
				if (value > rightmost_x) rightmost_x = value;
				else if (value < leftmost_x) leftmost_x = value;
				double value_y = cdt.triangle(*it)[i].y();
				if (value_y < bottom_y) bottom_y = value_y;
				else if (value_y > top_y) top_y = value_y;
			}
			counter++;
		}
		barys.close();
		for (auto it(triangles.begin()); it != triangles.end(); it++)
		{
			is_boundary[k] = false;
			if (CGAL::centroid(cdt.triangle(*it)).x() < 0.2) is_boundary[k] = true;
			for (int i = 0; i < 3; i++)
			{
				int index = k*3+i;
				auto opposite_tr = (*it)->neighbor(i);
				auto neighbor = find(triangles.begin(), triangles.end(), opposite_tr);
				if ((opposite_tr == nullptr) || (neighbor == triangles.end()))
				{
					neighbor_ids[index] = -1;
				}
				else
				{
					neighbor_ids[index] = distance(triangles.begin(), neighbor);
				}	
			}
			k++;
		}	
		k = 0;
		for (auto it(triangles.begin()); it != triangles.end(); it++)
		{
			Point barry = CGAL::centroid(cdt.triangle(*it));
			for (int i = 0; i < 3; i++)
			{
				Point corner1 = cdt.triangle(*it)[(i+1)%3];
				Point corner2 = cdt.triangle(*it)[(i+2)%3];
				Point center = Point((corner1.x() + corner2.x())/2.0, (corner1.y() + corner2.y())/2.0);
				double length = pow(pow(corner1.x()-corner2.x(), 2) + pow(corner1.y() - corner2.y(), 2), 0.5);
				double x1, y1;
				double x2, y2;
				if (fabs(corner1.y()-corner2.y())<=1e-6)
				{
					x1 = x2 = center.x();
					y1 = center.y() + length;
					y2 = center.y() - length;
				}
				else if (fabs(corner1.x()-corner2.x())<=1e-6)
				{
					y1 = y2 = center.y();
					x1 = center.x() + length;
					x2 = center.x() - length;
				}
				else
				{
					double n_x = (corner1.y()-corner2.y())/(corner2.x()-corner1.x());
					y1 = center.y() + pow(length*length/(1.0 + n_x*n_x),0.5);
					y2 = center.y() - pow(length*length/(1.0 + n_x*n_x),0.5);
					x1 = center.x() + n_x * (y1 - center.y());
					x2 = center.x() + n_x * (y2 - center.y());
				}
				Vector dist1 = Vector(x1 - barry.x(), y1 - barry.y());
				Vector dist2 = Vector(x2 - barry.x(), y2 - barry.y()); 
				Point dest = (dist2.squared_length() > dist1.squared_length()) ? Point(x2, y2) : Point(x1, y1);
				face_normals[2*(k*3+i)] = dest.x()-center.x();
				face_normals[2*(k*3+i)+1] = dest.y()-center.y();
			}
			k++;
		}
		k = 0;
		for (auto it(triangles.begin()); it != triangles.end(); it++)
		{
			tr_areas[k] = cdt.triangle(*it).area();
			k++;
		}
		//std::cout << "leaving the calc_dynamic_arrays " << std::endl;
	}

	void advection_2d::solver_small(double *&u, double const dt, double *&velocities, bool if_y, bool if_h, double t)
	{
		if (solver_type == 0)
		{
			double * interpolated_velocities;
			interpolated_velocities = new double [size * 3 * 2];
			solver_small_cpu(u, dt, velocities, interpolated_velocities, if_y, if_h, t);
			delete [] interpolated_velocities;
		}
		else
		{
			call_solver_small_cuda(u, dt, velocities, is_boundary, face_normals, neighbor_ids, tr_areas, xbarys, ybarys, size, if_y, if_h, t);
		}

	}

//mu_1 = 1.0, mu_2 = 1.0, epsilon_1 = 1.0, epsilon_2 = 2.25
	void advection_2d::solver_small_cpu(double *&u, double const &dt, double *&velocities, double *&interpolated_velocities, bool if_y, bool if_h, double t)
	{
		#pragma omp parallel for
		for (int j = 0; j < size; j++)
		{	
			//if(t==0.0) {int id = omp_get_thread_num();

        		//printf("thread: %i\n",id);}
			for (int k = 0; k < 3; k++)
			{
				if (neighbor_ids[j*3+k] != -1)
				{
					if (if_h && if_y)
					{
						if (is_boundary[j]!=is_boundary[neighbor_ids[j*3+k]] && ybarys[j] < 0.6 && ybarys[j] > 0.1 && ybarys[neighbor_ids[j*3+k]] < 0.6 && ybarys[neighbor_ids[j*3+k]] > 0.1)
						{
							interpolated_velocities[2*(j*3+k)] = (is_boundary[j] ? 1.0 : -1.0)*cos(25.0*(xbarys[neighbor_ids[j*3+k]] - 1.5*t))+velocities[2*j]+velocities[2*neighbor_ids[j*3+k]];
							interpolated_velocities[2*(j*3+k)+1] = velocities[2*j+1]+velocities[2*neighbor_ids[j*3+k]+1];
						}
						else
						{
							interpolated_velocities[2*(j*3+k)] = velocities[2*neighbor_ids[j*3+k]] + velocities[2*j];
							interpolated_velocities[2*(j*3+k)+1] = velocities[2*neighbor_ids[j*3+k]+1] + velocities[2*j+1];
						}
					}
					else
					{
						interpolated_velocities[2*(j*3+k)] = velocities[2*neighbor_ids[j*3+k]] + velocities[2*j];
						interpolated_velocities[2*(j*3+k)+1] = velocities[2*neighbor_ids[j*3+k]+1] + velocities[2*j+1];
					}
				}
				else
				{
					interpolated_velocities[2*(j*3+k)] = velocities[2*j]+0.0;
					interpolated_velocities[2*(j*3+k)+1] = velocities[2*j+1]+0.0;
				}
			}
		
		}
		#pragma omp parallel for
		for (int j = 0; j < size; j++)
		{
			//if (t==0.0) {int id = omp_get_thread_num();

        		//printf("thread: %i\n", id);}
			double temp = 0.0;
			for (int k = 0; k < 3; k++)
			{
				temp += interpolated_velocities[2*(3*j+k)] * face_normals[2*(j*3+k)] + interpolated_velocities[2*(3*j+k)+1] * face_normals[2*(j*3+k)+1]; 
			}
			if (!if_h) u[j] = u[j] - dt * (temp / tr_areas[j] + (if_y ? 0.5*pow(2.25, 0.5)*sigma(ybarys[j], 0.01, 7, 70)*u[j] : 0.5*pow(2.25, 0.5)*sigma(xbarys[j], 0.01, 15, 150)*u[j]) )/(1.0+0.5*dt*( if_y ? pow(2.25, 0.5)*sigma(ybarys[j], 0.01, 7, 70) : pow(2.25, 0.5)*sigma(xbarys[j], 0.01, 15, 150) ) );
			else u[j] = u[j] - dt * (temp / tr_areas[j] + 0.5*pow(1.0/2.25, 0.5)*sigma(ybarys[j], 0.01, 7, 70)*u[j] + 0.5*pow(1.0/2.25, 0.5)*sigma(xbarys[j], 0.01, 15, 150)*u[j])/(1.0+0.5*dt*( pow(1.0/2.25, 0.5)*sigma(ybarys[j], 0.01, 7, 70) + pow(1.0/2.25, 0.5)*sigma(xbarys[j], 0.01, 15, 150)));
		}
	}

	void advection_2d_mpi::solver_small(double *&u, double const &dt, double *&velocities, bool if_y, bool if_h, double t)
	{
		//cout << "entered solver, rank " << rank << endl; 
		vector<double> * send_data = new vector<double> [size_of_group];

		double ** receive_data = new double * [size_of_group];
		for (int i = 0; i < size_of_group; i++) receive_data[i] = new double [interfaces_size[i]*3];

		for (int k = 0; k < size; k++)
		{
			for (int i = 0; i < size_of_group; i++)
			{
				if (send_to[i][k])
				{
					send_data[i].push_back(u[k]);
					send_data[i].push_back(velocities[2*k]);
					send_data[i].push_back(velocities[2*k+1]);
				}
			}
		}
		//cout << "entered some other prt of solver, rank " << rank << endl;
		//MPI_Barrier(MPI_COMM_WORLD);
		
		MPI_Request request[size_of_group];
		//U EXCHANGE HERE
		for (int i = 0; i < size_of_group; i++)
		{
			//cout << i << " " << interfaces_size[i] << " " << send_data[i].size() << " " << rank << endl;
			double * pointer = send_data[i].data();
			if (interfaces_size[i]!=0) MPI_Isend(pointer, send_data[i].size(), MPI_DOUBLE, i, rank, MPI_COMM_WORLD, &request[i]);
		}

		for (int i = 0; i < size_of_group; i++)
		{
			if (interfaces_size[i]!=0) MPI_Recv(receive_data[i], interfaces_size[i]*3, MPI_DOUBLE, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		//MPI_Status status;
		for (int i = 0; i < size_of_group; i++)
		{
			if (interfaces_size[i]!=0) MPI_Wait(&request[i], MPI_STATUS_IGNORE);
		}

		//U EXCHANGE HERE END
		for (int i = 0; i < size_of_group; i++)
		{
			if (interfaces_size[i]!=0)
			{
				for (int k = 0; k < interfaces_size[i]; k++)
				{
					u[received_ids[i][k]] = receive_data[i][3*k];
					velocities[2*received_ids[i][k]] = receive_data[i][3*k+1];
					velocities[2*received_ids[i][k]+1] = receive_data[i][3*k+2];
				}
			}
		}
		for (int i = 0; i < size_of_group; i++) delete [] receive_data [i];
		delete [] receive_data;
		advection_2d::solver_small(u, dt, velocities, if_y, if_h, t);
		//cout << "finished solver, rank " << rank << endl;
		
	}

	void advection_2d_mpi::calc_mpi_arrays(CDT& cdt, vector<Face_handle>& local_triangles)
	{
		send_to = new bool * [size_of_group];
		for (int i = 0; i < size_of_group; i++) send_to[i] = new bool [size];
		for (int i = 0; i < size_of_group; i++) for (int j = 0; j < size; j++)
		{
			send_to[i][j] = false;
		}
		global_ids_reference = new int [size];

		int * counter = new int [size_of_group];
		for (int i =0; i < size_of_group; i++) counter[i] = 0;

		for (int j = 0; j < local_triangles.size(); j++)
		{
			Face_handle global_tr = get_find(cdt, local_triangles[j]);
			global_ids_reference[j] = get_distance(cdt, global_tr);	
			for (int num = 0; num< size_of_group; num++)
			{
				if (rank != num)
				{
					bool found = false;
					for (int k = 0; k < interfaces_size[num]; k++)
					{
						if (ids_to_send[num][k] == global_ids_reference[j])
						{
							found = true;
							break;
						}
					}
					if (found)
					{
						counter[num] ++;
						send_to[num][j] = true;
					}
				}
			}
		}
		//cout << "interface " << interfaces_size[0] << " " << interfaces_size[1] << " " << interfaces_size[2] << " " << interfaces_size[3] << " " << interfaces_size[4] << " "<< rank << endl;
		//cout << "result "<<  counter[0] << " " << counter[1] << " " << counter[2] << " " << counter[3] << " " << counter[4] <<" " << rank << " " << endl;
		received_ids = new int * [size_of_group];
		for (int i =0; i < size_of_group; i++) received_ids[i] = new int [interfaces_size[i]];

		MPI_Request request[size_of_group];
		MPI_Status status[size_of_group];
		for (int i = 0; i < size_of_group; i++)
		{
			if (i!=rank) MPI_Isend(ids_to_send[i], interfaces_size[i], MPI_INT, i, rank, MPI_COMM_WORLD, &request[i]);
		}
		for (int i = 0; i < size_of_group; i++)
		{
			if (i != rank)
			{
				int * buffer = new int [interfaces_size[i]];
				MPI_Recv(buffer, interfaces_size[i], MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
				for (int j = 0; j < size; j++)
				{
					for (int k = 0; k < interfaces_size[i]; k++)
					{
						if(global_ids_reference[j] == buffer[k]) received_ids[i][k] = j;
					}	
				}
				delete [] buffer;
			}
		}
	}
}
