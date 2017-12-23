# EcoSLIM
EcoSLIM private repo for development.  Will be public repo under /ParFlow when published


 SLIM2 is a Lagrangian, particle-tracking model for simulating subsurface and surface transport of reactive (such as microbial agents and metals) and non-reactive contaminants, diagnosing travel times, paths etc., which integrates
! seamlessly with ParFlow.
!
! Developed by: Reed Maxwell-August 2016 (rmaxwell@mines.edu)
!
! Contributors: Mohammad Danesh-Yazdi (danesh@mines.edu)
!
! Explore https://inside.mines.edu/~rmaxwell/maxwell_software.shtml
! for examples of ParFlow and SLIM (older version of SLIM2) working
! together.
!

! MASTER FORTRAN CODE

! slim2_main.f90: The main fortran code performing particke tracking
!                computations. To compile and execute on Mac OX
!		 operating ssytems:
!
!		 (i) first type the following in terminal:
!	 	     gfortran slim2_main.f90 -o compiler_name
!
!		 (ii) then typein terminal: ./compiler_name			
!
!--------------------------------------------------------------------
! INPUTS
!--------------------------------------------------------------------
! slimin.txt: Includes the domain's geometric information,
!             ParFlow timesteps and their total number, and particles
!             initial locations.
!
!--------------------------------------------------------------------
! SUBROUTINES
!--------------------------------------------------------------------
! pfb_read(arg1,...,arg5).f90: Reads a ParFlow .pfb output file and
!                              stores it in a matrix. Arguments
!                              in order are:
!
!                              - arg1: Name of the matrix in which
!                                      ParFlow .pfb is stored
!                              - arg2: Corresponding .pfb file name,
!                                      e.g., test.out.press.00100.pfb
!                              - arg3: Number of cells in x-direction
!                              - arg4: Number of cells in y-direction
!                              - arg5: Number of cells in z-direction
!
!--------------------------------------------------------------------
! OUTPUTS
!--------------------------------------------------------------------
! XXXX_log.txt:  Reports the domain's geometric information,
!                ParFlow's timesteps and their total number,
!                and particles initial condition. XXXX is the name of
!                the SLIM2 run already set in slimin.txt
!
! XXXX_particle.3D: Contains particles' trajectory information in
!                   space (i.e., X, Y, Z) and time (i.e., residence
!                   time). XXXX is the name of the SLIM2 run
!                   already set in slimin.txt
!
! XXXX_endparticle.txt: Contains the final X, Y, Z location of all
!                       particles as well as their travel time.
!                       XXXX is the name of the SLIM2 run already set
!                       in slimin.txt
!--------------------------------------------------------------------
