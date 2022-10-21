
#CXX := /share/home/wangw/gnu-need/mpich/bin/mpic++
CXX := mpiicpc
CC := mpiicc
#CC := /share/home/wangw/gnu-need/mpich/bin/mpicc
NETCDF := /share/home/wangw/intel-need/netcdfc
HDF5 := /share/home/wangw/intel-need/hdf5
ZLIB := /share/home/wangw/intel-need/zlib
SZ := /share/home/wangw/compression_comparison/install/SZ3

CFLAGS := -I../lib -I$(SZ)/include -I$(NETCDF)/include -I$(HDF5)/include -I$(ZLIB)/include  $(CFLAGS)

CXXFLAGS := -fast -std=c++0x $(CFLAGS)
#  -std=gnu++11

LDFLAGS := -L$(NETCDF)/lib -lnetcdf -L$(HDF5)/lib -lhdf5 -L$(SZ)/lib64 -lzstd -lstdc++ 

default: sz_mpi

all: sz_mpi

sz_mpi: \
	cJSON.o sz_mpi.o 
	$(CXX) -o $@ $^ $(LDFLAGS)

cJSON.o: ../lib/cJSON.c 
	${CC} -c -o $@ $(CXXFLAGS) $<

sz_mpi.o: sz_mpi.cpp 
	${CXX} -c -o $@ $(CXXFLAGS) $<

cleanall: 
	rm -f *.o sz_mpi
