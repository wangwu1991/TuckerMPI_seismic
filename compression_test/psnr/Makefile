
CXX := mpiicpc
CC := mpiicc
NETCDF := /share/home/wangw/intel-need/netcdfc
HDF5 := /share/home/wangw/intel-need/hdf5
ZLIB := /share/home/wangw/intel-need/zlib

CFLAGS := -I../lib -I$(NETCDF)/include -I$(HDF5)/include -I$(ZLIB)/include  $(CFLAGS)

CFLAGS := -fast -std=c99 $(CFLAGS)

LDFLAGS := -L$(NETCDF)/lib -lnetcdf -L$(HDF5)/lib -lhdf5  -lm 

default: psnr_mpi

all: psnr_mpi

psnr_mpi: \
	cJSON.o psnr_mpi.o 
	$(CC) -o $@ $^ $(LDFLAGS)
cJSON.o: ../lib/cJSON.c 
	${CC} -c -o $@ $(CFLAGS) $<
psnr_mpi.o: psnr_mpi.c 
	${CC} -c -o $@ $(CFLAGS) $<

cleanall: 
	rm -f *.o psnr_mpi
