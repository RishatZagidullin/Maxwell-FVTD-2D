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

namespace advection
{
	void initialize_arrays(int solver_type, int size, bool *& is_boundary, double *& face_normals, int *& neighbor_ids, double *& tr_areas, double *& xbarys, double *& ybarys);
	void dealloc_arrays(int solver_type, bool *& is_boundary, double *& face_normals, int *& neighbor_ids, double *& tr_areas, double *& xbarys, double *& ybarys);
	void call_solver_small_cuda(double *&u, double const dt, double *&velocities, bool * is_boundary, double * face_normals, int * neighbor_ids, double * tr_areas, double * xbarys, double * ybarys, int size, bool if_y, bool if_h, double t);
}
