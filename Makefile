NVPTH=/usr/local/nvidia/hpc_sdk/Linux_x86_64/23.7/compilers
export PATH := $(PATH):$(NVPTH)/bin
#LD_LIBRARY_PATH=$NVPTH/lib:/usr/local/nv_netcdf/lib64:$LD_LIBRARY_PATH
CC = nvcc
FC = nvfortran
F90 = nvfortran
OBJ_CORE = general.o netcdflib.o Definitions.o CPUcodes.o gpucodes.o DBscan.o Init.o Cells.o  mainsubs.o trj_analysis.o  ex-scan.o 

#FCOPTS = $(INC) -cuda -mavx2 -mno-avx512f -O2 -gpu=cc60,cc70,cc80 
FCOPTS = $(INC) -cuda -O2 -gpu=cc60,cc70,cc80 
F90OPTS =  -cuda -O2 -gpu=cc60,cc70,cc80
LKOPTS = -L/usr/local/nv_netcdf/lib64 -cuda -c++libs -gpu=cc80,cc70,cc60
LKLIBS = -lnetcdff
INC = -I /usr/local/nv_netcdf/include
#LKOPTS = -cuda -gpu=cc80,cc70,cc60


%.o : %.mod

.SUFFIXES : .f90

%.o: %.f90
	$(FC)  -c $(FCOPTS) $< -o $@

all: $(OBJ_CORE) 
	$(FC) $(LKOPTS) -o trj_analysis.exe $(OBJ_CORE)  $(LKLIBS)

ex-scan.o:
#	$(CC) -O4 --std c++14 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_80,code=sm_80 -c ex-scan.cu
	$(CC) -O4 --std c++14 -c ex-scan.cu

clean:
	/bin/rm -f *.o *.mod *.exe

