#include "example_meshes.h"

namespace example_meshes
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
		std::cout << "IT'S NOT FOUND" << std::endl;
		return count;
	}
	
	void mesh_mpi(CDT& cdt, double dr, double x_length, double y_length, int n_processes, double const * proportions, vector<Face_handle>& output, int ** & ids_to_send, int * & ids_size, int rank)
	{
		std::ofstream colours;
		if (rank==0) colours.open("colours.txt");

		Vertex_handle va = cdt.insert(Point(0.0,0.0));
  		Vertex_handle vb = cdt.insert(Point(x_length,0));
  		Vertex_handle vc = cdt.insert(Point(x_length, y_length));
  		Vertex_handle vd = cdt.insert(Point(0.0, y_length));
  		cdt.insert_constraint(va, vb);
  		cdt.insert_constraint(vb, vc);
  		cdt.insert_constraint(vc, vd);
  		cdt.insert_constraint(vd, va);
		CGAL::refine_Delaunay_mesh_2(cdt, Criteria(dr, 0.008));
		if (rank==0)
		{
			std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
			std::cout<<"TOTAL NUMBER OF FACES: " << cdt.number_of_faces()<<std::endl;
		}
		int * proc_labels = new int [cdt.number_of_faces()];
		for (int i = 0; i < cdt.number_of_faces(); i++) proc_labels[i] = 0;
		
		int i = 0;
		for (auto itt(cdt.finite_faces_begin()); itt != cdt.finite_faces_end(); itt++)
		{
			auto val = CGAL::centroid(cdt.triangle(itt));
			for (int j = 0; j < n_processes; j++)
			{
				if (val.x()>=proportions[j*4] && val.x()<proportions[j*4+1] && val.y()>=proportions[j*4+2] && val.y()<proportions[j*4+3]) proc_labels[i] = j+1;
			}
			i++;
		}
		cout << "finished one part" << endl;
		bool notstop = true;
		int counter = 0;
		while (notstop && counter < 5)
		{
			int k = -1;
			for (auto itt(cdt.finite_faces_begin()); itt != cdt.finite_faces_end(); itt++)
			{
				k++;
				//int c0 = get_distance(cdt, itt);
				//if (k!= c0) cout << "what the fuck" << endl;
				if ( (abs(proportions[4*(proc_labels[k]-1)]-CGAL::centroid(cdt.triangle(itt)).x()) < 0.01) || (abs(proportions[4*(proc_labels[k]-1)+1]-CGAL::centroid(cdt.triangle(itt)).x()) < 0.01) || (abs(proportions[4*(proc_labels[k]-1)+2]-CGAL::centroid(cdt.triangle(itt)).y()) < 0.01) || (abs(proportions[4*(proc_labels[k]-1)+3]-CGAL::centroid(cdt.triangle(itt)).y()) < 0.01) )
				{
				//cout << "this is how often it happens" << endl;
				cout << k << endl;
				if (proc_labels[k] > 0)
				{
					auto n1 = itt->neighbor(0);
					auto n2 = itt->neighbor(1);
					auto n3 = itt->neighbor(2);
					
					if (cdt.is_infinite(n1) || cdt.is_infinite(n2) || cdt.is_infinite(n3)) continue;

					int c1 = get_distance(cdt, n1);
					int c2 = get_distance(cdt, n2);
					int c3 = get_distance(cdt, n3);

					if (proc_labels[c1] == proc_labels[k] && proc_labels[c1] == proc_labels[c2]) ;
					else if (proc_labels[c2] == proc_labels[k] && proc_labels[c2] == proc_labels[c3]) ;
					else if (proc_labels[c3] == proc_labels[k] && proc_labels[c1] == proc_labels[c3]) ;
					else if (proc_labels[c1] == proc_labels[k] && proc_labels[c2] == proc_labels[c3]) 
					{
						 proc_labels[k] = proc_labels[c2];
					}
					else if (proc_labels[c2] == proc_labels[k] && proc_labels[c1] == proc_labels[c3])
					{
						 proc_labels[k] = proc_labels[c3];
					}
					else if (proc_labels[c3] == proc_labels[k] && proc_labels[c2] == proc_labels[c1])
					{
						proc_labels[k] = proc_labels[c1];
					}
				}
				else
				{
					std::cout << "THIS IS BAD";
				}
				}
				else
				{
					//cout << "this is how often it doesn't happen" << endl;
				}

			}
			notstop = false;
			counter++;
			std::cout << "second counter " << counter << std::endl;
		}
		/*int count1 = 0;
		int count2 = 0;
		int count3 = 0;
		int count4 = 0;
		int count5 = 0;
		double * p_new = new double [cdt.number_of_faces()];
		for (int k = 0; k < cdt.number_of_faces(); k++)
		{
			//p_new[k] = p_n[k];
			p_new[k] = ((float) proc_labels[k]) / ((float) n_processes);
		}
		for (int k = 0; k < cdt.number_of_faces(); k++)
		{
			int index = k;
			if ((p_new[index]) <= 0.25) colours << 0.0 << " " << ((p_new[index])/0.25)  << " " << 1.0 << " " << 1.0 << " "; 
			else if ((p_new[index]) <= 0.5) colours << 0.0 << " " << 1.0 << " " << ((0.5-(p_new[index]))/0.25) << " " << 1.0 << " ";
			else if ((p_new[index]) <= 0.75) colours << (((p_new[index])-0.5)/0.25) << " " << 1.0 << " " << 0.0 << " " << 1.0 << " ";
			else colours << 1.0 << " " << ((1.0-p_new[index])/0.25) << " " << 0.0 << " " << 1.0 << " ";	
			colours << endl;
		}
		delete [] p_new;
		cout << "finished drawing" <<endl;

			//if (rank ==0) colours << 1.0 << " " << 1.0 - ((float) proc_labels[i]) / ((float) n_processes) << " " << 1.0 - ((float) proc_labels[i]) / ((float) n_processes) << " " << 1.0 << std::endl;
			if (proc_labels[i]==0)
			{
				if(rank==0) colours << 1.0 << " " << 0.0 << " " << 1.0 << " " << 1.0 << std::endl;
			}
			if (proc_labels[i]==1)
			{
				if(rank==0) colours << 1.0 << " " << 0.0 << " " << 0.0 << " " << 1.0 << std::endl;
				count1++;
			}
			if (proc_labels[i]==2)
			{
				if(rank==0) colours << 0.0 << " " << 1.0 << " " << 0.0 << " " << 1.0 << std::endl;
				count2++;
			}
			if (proc_labels[i]==3)
			{
				if(rank==0) colours << 0.0 << " " << 0.0 << " " << 1.0 << " " << 1.0 << std::endl;
				count3++;
			}
			if (proc_labels[i]==4)
			{
				if(rank==0) colours << 0.0 << " " << 1.0 << " " << 1.0 << " " << 1.0 << std::endl;
				count4++;
			}
			if (proc_labels[i]==5)
			{
				if(rank==0) colours << 1.0 << " " << 1.0 << " " << 0.0 << " " << 1.0 << std::endl;
				count5++;
			}*/
		//std::cout << "count1: " << count1 << "count2: " << count2 << "count3: " << count3 << "count4: " << count4 << "count5: " << count5 << std::endl;
		colours.close();


		vector<int> temp_ids_to_send[n_processes][n_processes];
		vector<Face_handle> out [n_processes];


		int k = 0;
		for (auto itt(cdt.finite_faces_begin()); itt != cdt.finite_faces_end(); itt++)
		{
			if ( (abs(proportions[4*(proc_labels[k]-1)]-CGAL::centroid(cdt.triangle(itt)).x()) < 0.001) || (abs(proportions[4*(proc_labels[k]-1)+1]-CGAL::centroid(cdt.triangle(itt)).x()) < 0.001) || (abs(proportions[4*(proc_labels[k]-1)+2]-CGAL::centroid(cdt.triangle(itt)).y()) < 0.001) || (abs(proportions[4*(proc_labels[k]-1)+3]-CGAL::centroid(cdt.triangle(itt)).y()) < 0.001) )
			{
			cout << k << endl;
			//int dist = get_distance(cdt, itt);
			for (int num = 0; num < 3; num++) if (!cdt.is_infinite(itt->neighbor(num)) )
			{
				int val = proc_labels[get_distance(cdt, itt->neighbor(num))];
				if (proc_labels[k] != val ) temp_ids_to_send[proc_labels[k]-1][val-1].push_back(k);
			}
			}
			k++;
		}
		k = 0;
		for (auto itt(cdt.finite_faces_begin()); itt != cdt.finite_faces_end(); itt++)
		{
			for (int num = 0; num < n_processes; num++)
			{
				auto global_tr = find(temp_ids_to_send[num][rank].begin(), temp_ids_to_send[num][rank].end(), k);
				if (global_tr!=temp_ids_to_send[num][rank].end()) output.push_back(itt);
			}
			if (rank==proc_labels[k]-1) output.push_back(itt);
			cout << k++ << endl;
		}
		for (int i = 0; i < n_processes; i++)
		{
			ids_to_send[i] = new int [temp_ids_to_send[rank][i].size()];
			for (int j = 0; j < temp_ids_to_send[rank][i].size(); j++)
			{
				ids_to_send[i][j] = temp_ids_to_send[rank][i][j];
			}

			ids_size[i] = temp_ids_to_send[rank][i].size();
		}


	}
}


