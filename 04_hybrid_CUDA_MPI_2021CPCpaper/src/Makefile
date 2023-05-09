CUDA_INSTALL_PATH = /usr/local/cuda
MPI_INSTALL_PATH = /home/enuf_proj/mpich-install

NVCC = $(CUDA_INSTALL_PATH)/bin/nvcc -arch=compute_12 -code=compute_12 -DSINGLE_PRECISION
MPICC = $(MPI_INSTALL_PATH)/bin/mpicc -DSINGLE_PRECISION

LDFLAGS = -L$(MPI_INSTALL_PATH)/lib -L$(CUDA_INSTALL_PATH)/lib64 
LIB = -lcuda -lcudart -lcufft -lm -lmpich -lpthread

CUFILES = cunfft.cu
CFILES = util.c hpnfft_forward.c hpnfft_backward.c TEST.c
OBJECTS = cunfft.o util.o hpnfft_forward.o hpnfft_backward.o TEST.o
EXECNAME = CUDA_MPI_TEST

all:
	$(NVCC) -c $(CUFILES)
	$(MPICC) -c $(CFILES)
	$(MPICC) -o $(EXECNAME) $(OBJECTS) $(LIB) $(LDFLAGS)

clean:	
	rm -f $(OBJECTS)

