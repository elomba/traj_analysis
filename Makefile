NVPTH=/usr/local/nvidia/hpc_sdk/Linux_x86_64/23.7/compilers
export PATH := $(PATH):$(NVPTH)/bin
CC = nvcc
FC = nvfortran
F90 = nvfortran
OBJ_CORE = general.o netcdflib.o Definitions.o CPUcodes.o gpucodes.o DBscan.o Init.o Cells.o  mainsubs.o trj_analysis.o  ex-scan.o 

# Remove avx512 instructions to avoid 'Illegal instruction' errors in older processors
FCOPTS = $(INC) -cuda -fast -mavx2 -mno-avx512f -gpu=cc60,cc70,cc80 
LKOPTS = -L/usr/local/nv_netcdf/lib64 -cuda -c++libs -gpu=cc80,cc70,cc60
LKLIBS = -lnetcdff
INC = -I /usr/local/nv_netcdf/include


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

