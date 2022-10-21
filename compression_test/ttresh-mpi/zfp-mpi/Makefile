
# CXX := mpiicpc
CC := mpiicc
NETCDF := /share/home/wangw/intel-need/netcdfc
HDF5 := /share/home/wangw/intel-need/hdf5
ZLIB := /share/home/wangw/intel-need/zlib
ZFP := /share/home/wangw/compression_comparison/install/zfp

CFLAGS := -I../lib -I$(ZFP)/include -I$(NETCDF)/include -I$(HDF5)/include -I$(ZLIB)/include  $(CFLAGS)

CXXFLAGS := -fast -std=c99  $(CFLAGS)
# -std=c++0x -std=gnu++11

LDFLAGS := -L$(NETCDF)/lib -lnetcdf -L$(HDF5)/lib -hdf5 -L$(ZFP)/lib  -lm -lzfp 

default: zfp_mpi

all: zfp_mpi

zfp_mpi: \
	cJSON.o zfp_mpi.o 
	$(CC) -o $@ $^ $(LDFLAGS)

cJSON.o: ../lib/cJSON.c 
	${CC} -c -o $@ $(CXXFLAGS) $<

zfp_mpi.o: zfp_mpi.c
	${CC} -c -o $@ $(CXXFLAGS) $<

cleanall: 
	rm -f *.o zfp_mpi
