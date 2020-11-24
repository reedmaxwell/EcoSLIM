This is the multi-GPU version.

It can be built by:
pgf90 -o EcoSLIM.exe -mp ecoslim_cuda_md-test.cuf *.f90 EcoSLIM_CUDA-omp.cuf

If you use multi-GPU, please set:
export OMP_NUM_THREADS=n
where n is the number of GPU.

If you use only one GPU, please also set:
export OMP_NUM_THREADS=1

An example of CUDA environment:
Ubuntu 18.04.5
PGI 19.10
CUDA 10.2

An example of slimin.txt is in the Example folder
This branch has the same capability with the master branch.
One parameter is different from the master branch:
line 15 in slimin.txt: it should be the current timestep minus 1.
For example, if it is the first run of your simulation and you start from time=1.0, you should set it as 0


If there are any problems or you need technical support, please contact: chen_yang@princeton.edu
