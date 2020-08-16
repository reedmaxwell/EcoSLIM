objects = \
	ecoslim_cuda_m.o  \
	ran1.mod.o \
	pfb_read.o \
	vtk_write.o \
	vtk_write_points.o


.SUFFIXES:
.SUFFIXES:	.mod .o .f90

EcoSLIM.exe : $(objects)
		pgf90 -o EcoSLIM.exe $(objects) EcoSLIM_CUDA.cuf
		
.f90.o :
		pgf90 -c $<
		
ecoslim_cuda_m.o : ecoslim_cuda_m.cuf
		pgf90 -c ecoslim_cuda_m.cuf		

clean :
		rm -f EcoSLIM.exe $(objects) EcoSLIM_CUDA.o *.mod

.PHONY : debug
	debug :
		pgf90 -gDD -ffpe-trap=zero,overflow,underflow -fbacktrace -fbounds-check -o EcoSLIM.exe $(objects)
