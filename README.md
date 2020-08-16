# EcoSLIM
This is the MPI version.
Please make the Makefile to build the EcoSLIM.
Run example:
mpirun -np 20 ./EcoSLIM.exe
20 is the number of processes used.
If there is error, please try:
mpirun -np 20 ./EcoSLIM.exe --mca mpi_cuda_support 0

