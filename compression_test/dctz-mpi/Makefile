
# CXX := mpiicpc
CC := mpiicc
NETCDF := /share/home/wangw/intel-need/netcdfc
HDF5 := /share/home/wangw/intel-need/hdf5
ZLIB := /share/home/wangw/intel-need/zlib
fftw := /share/home/wangw/intel-need/fftw3
# SZ := /media/ww/data/paper1_prepare/compression_comparison/SZ
# ZSTD := /media/ww/data/paper1_prepare/compression_comparison/SZ3-initial/zstd
# ZFP := /media/ww/data/paper1_prepare/compression_comparison/zfp

CFLAGS := -I../lib -I$(NETCDF)/include -I$(HDF5)/include -I$(ZLIB)/include -I$(fftw)/include  $(CFLAGS)

CXXFLAGS := -fast -std=c99 -Wall -Wno-unused -DUSE_TRUNCATE -DUSE_QTABLE $(CFLAGS)
# -std=c++0x -std=gnu++11

LDFLAGS := -L$(NETCDF)/lib -lnetcdf -L$(HDF5)/lib -hdf5 -L$(fftw)/lib  -lm -lfftw3 -lfftw3f -lz -lpthread  -DUSE_TRUNCATE -DUSE_QTABLE

default: dctz_mpi

all: dctz_mpi

dctz_mpi: \
	cJSON.o dct-float.o dct.o util.o binning.o dctz-decomp-lib.o dctz-comp-lib.o dctz_mpi.o 
	$(CC) -o $@ $^ $(LDFLAGS)

cJSON.o: ../lib/cJSON.c 
	${CC} -c -o $@ $(CXXFLAGS) $<
dct-float.o: dct-float.c 
	${CC} -c -o $@ $(CXXFLAGS) $<
dct.o: dct.c 
	${CC} -c -o $@ $(CXXFLAGS) $<
util.o: util.c 
	${CC} -c -o $@ $(CXXFLAGS) $<
binning.o: binning.c 
	${CC} -c -o $@ $(CXXFLAGS) $<
dctz-decomp-lib.o: dctz-decomp-lib.c 
	${CC} -c -o $@ $(CXXFLAGS) $<
dctz-comp-lib.o: dctz-comp-lib.c 
	${CC} -c -o $@ $(CXXFLAGS) $<
dctz_mpi.o: dctz-mpi.c  
	${CC} -c -o $@ $(CXXFLAGS) $< 

cleanall: 
	rm -f *.o dctz_mpi
