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

If there are any problems or you need technique support, please contact: chen_yang@princeton.edu
