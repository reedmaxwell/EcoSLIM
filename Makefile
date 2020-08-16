objects = \
	ran1.mod.o \
	EcoSLIM_MPI.o  \
	pfb_read.o \
	vtk_write.o \
	vtk_write_points.o


.SUFFIXES:
.SUFFIXES:  .mod .o .f90

EcoSLIM.exe : $(objects)
		mpif90 -O3 -o EcoSLIM.exe $(objects)

.f90.o :
		mpif90 -O3 -c $<

clean :
		rm -f EcoSLIM.exe $(objects) *.mod
		
.PHONY : debug
	debug :
		mpif90 -gDD -ffpe-trap=zero,overflow,underflow -fbacktrace -fbounds-check -o EcoSLIM.exe $(objects)
