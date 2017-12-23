objects = \
	EcoSLIM.o  \
	pfb_read.o


.SUFFIXES:
.SUFFIXES:  .o .f90 

EcoSLIM.exe : $(objects)
		gfortran -fopenmp  -O3 -o EcoSLIM.exe $(objects)

.f90.o :
		gfortran -fopenmp -O3 -c $<

.PHONY : clean
	clean : 
		rm -f EcoSLIM.exe $(objects)

debug :
		gfortran -gDD -ffpe-trap=zero,overflow,underflow -fbacktrace -fbounds-check -o EcoSLIM.exe $(objects)
 
