# TuckerMPI_seismic
parallel Tucker tensor decomposition used for seismic data compression

# Attension
the code is forked from https://gitlab.com/tensors/TuckerMPI

# Major Changes
the code has modified into a function that can be called by other program to achieve online compression, so it can not be used by independent compilation; if you want a independent version of TuckerMPI, you can refer to https://gitlab.com/tensors/TuckerMPI;

# Acknowledge
We would like to thank that Grey Ballard, Tamara Kolda, Hemanth Kolla, Woody Austin, Alicia Klinvex, Casey Battaglino.

# Some codes for compression test 
These code are for manuscript "Compression of Seismic Forward Modeling Wavefield Using TuckerMPI", including ZFP, SZ2, SZ3, DCTZ, and TTHRESH. The users need to install those compressors on the computer in advance, then modify the path in the "Makefile" and compile them with MPI compiler. 

# Data
The test data have been uploaded to https://drive.google.com/drive/folders/1N1L2802-9LSEMRUhXpjmCutpLC9DEAOk?usp=sharing




