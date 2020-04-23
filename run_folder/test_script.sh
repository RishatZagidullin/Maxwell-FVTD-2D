#!/bin/bash -l
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --gres=gpu:1
#SBATCH --time=0-0:10:0
#SBATCH -p gpu_devel

module load gpu/cuda-9.2
module load compilers/gcc-7.3.0
module load mpi/openmpi-3.1.4

echo $LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/trinity/home/r.zagidullin/CGAL/lib64
echo $LD_LIBRARY_PATH

mpirun 2D_MPI_CUDA.exe input_2d_measure
