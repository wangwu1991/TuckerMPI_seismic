
CXX := mpiicpc
CC := mpiicc
NETCDF := /share/home/wangw/intel-need/netcdfc
HDF5 := /share/home/wangw/intel-need/hdf5
ZLIB := /share/home/wangw/intel-need/zlib
TTHRESH := ../src
Eigen := ../external


CFLAGS := -I./lib -I$(TTHRESH) -I$(Eigen) -I$(NETCDF)/include -I$(HDF5)/include -I$(ZLIB)/include  $(CFLAGS)

CXXFLAGS := -fast -std=c++11  $(CFLAGS)
#  -std=gnu++11

LDFLAGS := -L$(NETCDF)/lib -lnetcdf -L$(HDF5)/lib -lhdf5  -lm 

default: ttresh_mpi

all: ttresh_mpi

ttresh_mpi: \
	cJSON.o ttresh_mpi.o 
	$(CXX) -o $@ $^ $(LDFLAGS)

cJSON.o: lib/cJSON.c 
	${CXX} -c -o $@ $(CXXFLAGS) $<

ttresh_mpi.o: ttresh_mpi.cpp 
	${CXX} -c -o $@ $(CXXFLAGS) $<

cleanall: 
	rm -f *.o ttresh_mpi
