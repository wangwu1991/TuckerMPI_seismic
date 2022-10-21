
CXX := mpiicpc
CC := mpiicc
NETCDF := /share/home/wangw/intel-need/netcdfc
HDF5 := /share/home/wangw/intel-need/hdf5
ZLIB := /share/home/wangw/intel-need/zlib
SZ := /share/home/wangw/compression_comparison/install/SZ2

CFLAGS := -I../lib -I$(SZ)/include -I$(NETCDF)/include -I$(HDF5)/include -I$(ZLIB)/include  $(CFLAGS)

CFLAGS := -fast -std=c99 $(CFLAGS)
# -std=c++0x -std=gnu++11 $(SZ)/lib/libSZ.a $(SZ)/lib/libzlib.a $(SZ)/lib/libzstd.a -Wl,-rpath,$(SZ)/lib 

LDFLAGS := -L$(NETCDF)/lib -lnetcdf -L$(HDF5)/lib -lhdf5 $(SZ)/lib/libSZ.a $(SZ)/lib/libzlib.a $(SZ)/lib/libzstd.a -Wl,-rpath,$(SZ)/lib  -lm 

default: sz_mpi

all: sz_mpi

sz_mpi: \
	cJSON.o sz_mpi.o 
	$(CC) -o $@ $^ $(LDFLAGS)

cJSON.o: ../lib/cJSON.c 
	${CC} -c -o $@ $(CFLAGS) $<

sz_mpi.o: sz_mpi.c 
	${CC} -c -o $@ $(CFLAGS) $<

cleanall: 
	rm -f *.o sz_mpi
