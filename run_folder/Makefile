MPICC = mpicxx
NVCC = nvcc

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/trinity/home/r.zagidullin/CGAL/lib64
#export CPATH=$CPATH:/trinity/home/r.zagidullin/CGAL/include

CUDA_INSTALL_PATH = /usr/local/cuda#MSI
#CUDA_INSTALL_PATH = ${CUDA_PATH}
CUDACFLAGS=-I${CUDA_INSTALL_PATH}/include
CUDALDFLAGS=-L${CUDA_INSTALL_PATH}/lib64
CUDA_LIB = -lcudart


MPI_INSTALL_PATH = /usr/lib/x86_64-linux-gnu/openmpi#MSI and LENOVO
#MPI_INSTALL_PATH = ${MPI_HOME}
MPICFLAGS=-I${MPI_INSTALL_PATH}/include
MPILDFLAGS=-L${MPI_INSTALL_PATH}/lib
MPI_LIB = -lmpi -lmpi_cxx



LIB  = -L../fade2d/fade_lib -lgomp -lm -lgmp -lCGAL

GCC_COPT = -fopenmp
NVCC_COPT = -Xcompiler -fopenmp

VPATH = ../advection:../

SOURCES = example_meshes.cpp transport.cpp main_2d_mpi.cpp
OBJECTS = $(SOURCES:.cpp=_cpp.o)

CUDA_SOURCES = transport.cu helper.cu
CUDA_OBJECTS = $(CUDA_SOURCES:.cu=_cuda.o)

all: 2D_MPI_CUDA.exe

2D_MPI_CUDA.exe: $(OBJECTS) $(CUDA_OBJECTS)
	$(NVCC) $^ -o $@ -L/trinity/home/r.zagidullin/CGAL/lib64 $(LIB) $(NVCC_COPT) $(MPILDFLAGS) $(MPI_LIB)

%_cpp.o: %.cpp
	$(MPICC) -c $< -o $@ $(GCC_COPT) $(CUDACFLAGS)

%_cuda.o: %.cu
	$(NVCC) -c $< -o $@ $(NVCC_COPT) $(MPICFLAGS)

clean:
	rm -rf *.o
	rm -rf *.exe
