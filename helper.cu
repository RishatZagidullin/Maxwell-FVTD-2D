#include "helper.cuh"


void initialize_cuda_memory(double * &E_x, double * &E_y, double * &H_z, double * &coefs_E_x, double * &coefs_E_y, double * &coefs_H_z, double * &mu, double * &epsilon, int size)
{
	cudaMallocManaged((void **) &E_y, size*sizeof(double));
	cudaMallocManaged((void **) &E_x, size*sizeof(double));
	cudaMallocManaged((void **) &H_z, size*sizeof(double));
	cudaMallocManaged((void **) &coefs_E_x, size*2*sizeof(double));
	cudaMallocManaged((void **) &coefs_E_y, size*2*sizeof(double));
	cudaMallocManaged((void **) &coefs_H_z, size*2*sizeof(double));
	cudaMallocManaged((void **) &mu, size*sizeof(double));
	cudaMallocManaged((void **) &epsilon, size*sizeof(double));
}

void set_device(int rank)
{
	cudaSetDevice(rank);
}

void call_update_coefs_cuda(double * H_z, double * E_x, double * E_y, double *coefs_H_z, double *coefs_E_x, double *coefs_E_y, double *mu, double *epsilon, double size)
{
	dim3 block(32);
	dim3 grid_big((size+block.x-1)/block.x);
	update_coefs_cuda<<<grid_big, block>>>(H_z, E_x, E_y, coefs_H_z, coefs_E_x, coefs_E_y, mu, epsilon, size);
	cudaDeviceSynchronize();
}

__global__ void update_coefs_cuda(double * H_z, double * E_x, double * E_y, double *coefs_H_z, double *coefs_E_x, double *coefs_E_y, double *mu, double *epsilon, double size)
{
	int numThreads = blockDim.x * blockDim.y * gridDim.x * gridDim.y;
	int global_id = (threadIdx.y + blockIdx.y * blockDim.y)*blockDim.x*gridDim.x + (threadIdx.x + blockIdx.x * blockDim.x);
	for (int p = global_id; p < size; p+=numThreads)
	{	
		coefs_H_z[2*p] = 1.0/mu[p]*E_y[p];
		coefs_H_z[2*p+1] = -1.0/mu[p]*E_x[p];
		coefs_E_x[2*p] = 0.0;
		coefs_E_y[2*p] = 1.0/epsilon[p]*H_z[p];
		coefs_E_x[2*p+1] = -1.0/epsilon[p]*H_z[p];
		coefs_E_y[2*p+1] = 0.0;
	}
}
