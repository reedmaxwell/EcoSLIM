!--------------------------------------------------------------------
! **EcoSLIM** is a Lagrangian, particle-tracking that simulates advective
! and diffusive movement of water parcels.  This code can be used to
! simulate age, diagnosing travel times, source water composition and
! flowpaths.  It integrates seamlessly with **ParFlow-CLM**.
!
! Developed by: Reed Maxwell-August 2016 (rmaxwell@mines.edu)
!
! Contributors: Laura Condon (lecondon@email.arizona.edu)
!               Mohammad Danesh-Yazdi (danesh@sharif.edu)
!               Lindsay Bearup (lbearup@usbr.gov)
!
! released under GNU LPGL, see LICENSE file for details
!
!--------------------------------------------------------------------
! MAIN FORTRAN CODE
!--------------------------------------------------------------------
! EcoSLIM_main.f90: The main fortran code performing particle tracking
!                 Use Makefile provided to build.
!
!
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
!
!--------------------------------------------------------------------
! CODE STRUCTURE
!--------------------------------------------------------------------
! (1) Define variables
!
! (2) Read inputs, set up domain, write the log file, and
!     initialize particles,
!
! (3) For each timestep, loop over all particles to find and
!     update their new locations
!--------------------------------------------------------------------


program EcoSLIM
use ran_mod
use mpi
implicit none
!--------------------------------------------------------------------
! (1) Define variables
!--------------------------------------------------------------------

real*8,allocatable::P(:,:)
        ! P = Particle array [np,attributes]
        ! np = Number of particles
        ! P(np,1) = X coordinate [L]
        ! P(np,2) = Y coordinate [L]
        ! P(np,3) = Z coordinate [L]
        ! P(np,4) = Particle residence time [T]
        ! P(np,5) = Saturated particle residence time [T]
        ! P(np,6) = Particle mass; assigned via preciptiation or snowmelt rate (Evap_Trans*density*volume*dT)
        ! P(np,7) = Particle source (1=IC, 2=rain, 3=snowmelt...)
        ! P(np,8) = Particle Status (1=active, 0=inactive)
        ! P(np,9) = concentration
        ! P(np,10) = Exit status (1=outflow, 2=ET...)

!@ RMM, why is this needed?
real*8,allocatable::PInLoc(:,:)
        ! PInLoc(np,1) = Particle initial X location
        ! PInLoc(np,2) = Particle initial Y location
        ! PInLoc(np,3) = Particle initial Z location

real*8,allocatable::Vx(:,:,:)
real*8,allocatable::Vy(:,:,:)
real*8,allocatable::Vz(:,:,:)
        ! Vx = Velocity x-direction [nx+1,ny,nz] -- ParFlow output
        ! Vy = Velocity y-direction [nx,ny+1,nz] -- ParFlow output
        ! Vz = Velocity z-direction [nx,ny,nz+1] -- ParFlow output

real*8,allocatable::C(:,:,:,:)
real*8,allocatable::C_p(:,:,:,:)
        ! Concentration array, in i,j,k with l (first index) as consituent or
        ! property.  These are set by user at runtime using input

!real*8,allocatable::ET_grid(:,:,:,:)
        ! ET array, in i,j,k with l (first index) as consituent or property

CHARACTER*20,allocatable:: conc_header(:)
        ! name for variables written in the C array above.  Dimensioned as l above.
real*8,allocatable::Time_Next(:)
        ! Vector of real times at which ParFlow dumps outputs

real*8,allocatable::dz(:), Zt(:)
        ! Delta Z values in the vertical direction
        ! Elevations in z-direction in local coordinates

real*8,allocatable::Sx(:,:)  ! Sx: Slopes in x-direction (not used)
real*8,allocatable::Sy(:,:)  ! Sy: Slopes in y-direction (not used)

real*8,allocatable::Saturation(:,:,:)    ! Saturation (read from ParFlow)
real*8,allocatable::Porosity(:,:,:)      ! Porosity (read from ParFlow)
real*8,allocatable::EvapTrans(:,:,:)     ! CLM EvapTrans (read from ParFlow, [1/T] units)
real*8,allocatable::CLMvars(:,:,:)     ! CLM Output (read from ParFlow, following single file
                                       ! CLM output as specified in the manual)
real*8, allocatable::Pnts(:,:), DEM(:,:) ! DEM and grid points for concentration output

integer Ploc(3)
        ! Particle's location whithin a cell

integer nx, nnx, ny, nny, nz, nnz, nztemp
        ! number of cells in the domain and cells+1 in x, y, and z directions

integer np_ic, np, np_active, np_active2, icwrite, jj, npnts, ncell, npout
        ! number of particles for intial pulse IC, total, and running active

integer nt, n_constituents
        ! number of timesteps ParFlow; numer of C vectors written for VTK output

real*8  pfdt, advdt(3)
        ! ParFlow timestep value, advection timestep for each direction
        ! for each individual particle step; used to chose optimal particle timestep

integer pft1, pft2, tout1, pfnt, n_cycle
        ! parflow start and stop file numbers number of ParFlow timesteps
        ! flag specifying the number for the first output write (0= start with pft1)
        ! number of timestep cycles

real*8  Time_first
        ! initial timestep for Parflow ((pft1-1)*pfdt)

integer kk
        ! Loop counter for the time steps (pfnt)
integer pfkk, outkk
       ! Counter for the file numbers starts at pft1
       ! Counter for the output writing

integer ii
        ! Loop counter for the number of particles (np)
integer iflux_p_res, i_added_particles
        ! Number of particles per cell for flux input and
        ! number of particle added total for precip input per
        ! timestep
integer i, j, k, l, ik, ji, m, ij, nzclm, nCLMsoil
integer itime_loc
        ! Local indices / counters
integer*4 ir

character*200 runname, filenum, filenumout, pname, fname, vtk_file, DEMname
        ! runname = SLIM runname
        ! filenum = ParFlow file number
        ! filenumout = File number for Ecoslim writing
        ! pname = ParFlow output runname
        ! fname = Full name of a ParFlow's output
        ! vtk_file = concentration file
        ! DEMname = DEM file name

real*8 Clocx, Clocy, Clocz, Z, maxz
        ! The fractional location of each particle within it's grid cell
        ! Particle Z location

real*8 V_mult
        ! Multiplier for forward/backward particle tracking
        ! If V_mult = 1, forward tracking
        ! If V_mult = -1, backward tracking

logical clmtrans, clmfile
        ! logical for mode of operation with CLM, will add particles with P-ET > 0
        ! will remove particles if ET > 0
        ! clmfile governs reading of the full CLM output, not just evaptrans

real*8 dtfrac
        ! fraction of dx/Vx (as well as dy/Vy and dz/Vz) assuring
        ! numerical stability while advecting a particle to a new
        ! location.

real*8 Xmin, Xmax, Ymin, Ymax, Zmin, Zmax
        ! Domain boundaries in local / grid coordinates. min values set to zero,
        ! DEM is read in later to output to Terrain Following Grid used by ParFlow.

real*8 dx, dy
        ! Domain's number of cells in x and y directions

real*8 Vpx, Vpy, Vpz
        ! Particle velocity in x, y, and z directions

real*8 particledt, delta_time
        ! The time it takes for a particle to displace from
        ! one location to another and the local particle from-to time
        ! for each PF timestep

real*8  mean_age, mean_comp, mean_mass, total_mass
        ! mean age and composition and mass of all particles in domain

real*8 local_flux, et_flux, water_vol, Zr, z1, z2, z3
        ! The local cell flux convergence
        ! The volumetric ET flux
        ! The availble water volume in a cell
        ! random variable

real*8 Xlow, Xhi, Ylow, Yhi, Zlow, Zhi
        ! Particles initial locations i.e., where they are injected
        ! into the domain.

! density of water (M/L3), molecular diffusion (L2/T), fractionation
real*8 denh2o, moldiff, Efract  !, ran1

! time history of ET, time (1,:) and mass for rain (2,:), snow (3,:),
! PET balance is water balance flux from PF accumulated over the domain at each
! timestep
real*8, allocatable::ET_age(:,:), ET_comp(:,:), PET_balance(:,:), ET_mass(:)
real*8, allocatable::Out_age(:,:), Out_mass(:,:), Out_comp(:,:)
integer, allocatable:: ET_np(:), Out_np(:)

real*8::ET_age_p, ET_comp_p(3), ET_mass_p
real*8::Out_age_p, Out_mass_p, Out_comp_p(3)
integer:: ET_np_p, Out_np_p

real*8  ET_dt, DR_Temp
! time interval for ET
! integer counters and operators.
! the first set are used for total run timing the latter for component timing
integer*8  Total_time1, Total_time2, t1, t2, IO_time_read, IO_time_write, parallel_time
integer*8  sort_time, T3, T4
!! integers for writing C or point based output
integer ipwrite, ibinpntswrite
! integers for writing gridded ET outputs
integer etwrite

!! IO control
!! ipwrite controls an ASCII, .3D particle file not recommended due to poor performance
!! this is left as a compiler option, currently disabled
!!!!!!!
!! icwrite controls VTK, binary grid based output where particle masses, concentrations,
!! ages are mapped to a grid and written every N timesteps.  This is the most effiecient output
!! but loses some accuracy or flexbility becuase individual particle locations are aggregated
!! to the grid
!!!!!!!!
!! ibinpntswrite controls VTK, binary output of particle locations and attributes.  This is much faster
!! than the .3D ASCII output but is still much slower than grid based output.  It provides the most
!! information as particle locations are preserved

integer rank, ierr, pp, status(MPI_STATUS_SIZE)
integer, allocatable::np_p(:) 
integer BCAST_int(16)
real*8  BCAST_float(5)

interface
  SUBROUTINE vtk_write(time,x,conc_header,ixlim,iylim,izlim,icycle,n_constituents,Pnts,vtk_file)
  real*8                 :: time
  REAL*8    :: x(:,:,:,:)
  CHARACTER (LEN=20)     :: conc_header(:)
  INTEGER*4 :: ixlim
  INTEGER*4 :: iylim
  INTEGER*4 :: izlim
  REAL*8                 :: dx
  REAL*8                 :: dy
  REAL*8                 :: dz(izlim)
  REAL*8                 :: Pnts(:,:)
  INTEGER                :: icycle
  INTEGER                :: n_constituents
  CHARACTER (LEN=200)    :: vtk_file
end subroutine vtk_write

SUBROUTINE vtk_write_points(P,np_active, np,icycle,vtk_file,dx,dy,nx,ny,maxZ,DEM)
REAL*8    :: P(:,:)
INTEGER                :: icycle
INTEGER*4              :: np_active
INTEGER*4              :: np
INTEGER                :: n_constituents
REAL*8                 :: dx
REAL*8                 :: dy
INTEGER*4              :: nx
INTEGER*4              :: ny
REAL*8                 :: maxZ
REAL*8                 :: DEM(:,:)
CHARACTER (LEN=200)    :: vtk_file
end subroutine vtk_write_points

end interface

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, pp, ierr)

allocate(np_p(pp))

if(rank == 0) then

!Set up timing
Total_time1 = 0
Total_time2 = 0
t1 = 0
t2 = 0
IO_time_read = 0
IO_time_write = 0
parallel_time = 0
sort_time = 0

        call system_clock(Total_time1)

!--------------------------------------------------------------------
! (2) Read inputs, set up domain, write the log file, and
! initialize particles
!--------------------------------------------------------------------

! Note: The following file numbers refer to
!
!       - #10: slimin.txt
!       - #11: runname_log.txt
!       - #12: runname_particle.3D (visualizes particles in VisIT)
!       - #13: runname_endparticle.txt

call system_clock(T1)

! open SLIM input .txt file
open (10,file='slimin.txt')

! read SLIM run name
read(10,*) runname

! read ParFlow run name
read(10,*) pname

! read DEM file name
read(10,*) DEMname

! open/create/write the output log.txt file. If doesn't exist, it's created.
open(11,file=trim(runname)//'_log.txt')
write(11,*) '### EcoSLIM Log File'
write(11,*)
write(11,*) 'run name:',trim(runname)
write(11,*)
write(11,*) 'ParFlow run name:',trim(pname)
write(11,*)
if (DEMname /= '') then
write(11,*) 'ParFlow DEM name:',trim(DEMname)
else
write(11,*) 'Not reading ParFlow DEM'
end if
write(11,*)

! read domain number of cells and number of particles to be injected
read(10,*) nx
read(10,*) ny
read(10,*) nz
! read in number of particles for IC (if np_ic = -1 then restart from a file)
read(10,*) np_ic

! read in the number of particles total
read(10,*) np

! check to make sure we don't assign more particles for IC than we have allocated
! in total
if (np_ic > np) then
write(11,*) ' warning NP_IC greater than IC'
np = np_ic
end if

! write nx, ny, nz, and np in the log file
write(11,*) ' Grid information'
write(11,*) 'nx:',nx
write(11,*) 'ny:',ny
write(11,*) 'nz:',nz
write(11,*)
write(11,*) ' Particle IC Information'
write(11,*) 'np IC:',np_ic
if (np_ic == -1) &
write(11,*)  'Reading particle restart file:',trim(runname)//'_particle_restart.bin'
write(11,*) 'np:',np


! allocate P, Sx, dz, Vx, Vy, Vz, Saturation, and Porosity arrays
allocate(P(np,10))
P(1:np,1:6) = 0.0    ! clear out all particle attributes
P(1:np,7:9) = 1.0  ! make all particles active to start with and original from 1 = GW/IC
P(1:np,10) = 0.0   ! no exit

nCLMsoil = 10 ! number of CLM soil layers over the root zone
nzclm = 13+nCLMsoil ! CLM output is 13+nCLMsoil layers for different variables not domain NZ,
           !  e.g. 23 for 10 soil layers (default) and 17 for 4 soil layers (Noah soil
           ! layer setup)

!  number of things written to C array, hard wired at 2 now for Testing
!n_constituents = 5
n_constituents = 9
!allocate arrays
allocate(PInLoc(np,3))
allocate(Sx(nx,ny),Sy(nx,ny), DEM(nx,ny))
allocate(dz(nz), Zt(0:nz))
allocate(CLMvars(nx,ny,nzclm))
allocate(C(n_constituents,nx,ny,nz))
!allocate(ET_grid(1,nx,ny,nz))
allocate(conc_header(n_constituents))
C = 0.0d0

! read dx, dy as scalars
read(10,*) dx
read(10,*) dy

! read dz as an array
read(10,*) dz(1:nz)

! read in (constant for now) ParFlow dt
read(10,*) pfdt

! read in parflow start and stop times
read(10,*) pft1
read(10,*) pft2
read(10,*) tout1
read(10,*) n_cycle

pfnt=n_cycle*(pft2-pft1+1)

if (tout1 == 0) then
   outkk = pft1
   tout1 = pft1
else
  outkk = tout1
end if

! set ET DT to ParFlow one and allocate ET arrays accordingly
ET_dt = pfdt
allocate(ET_age(pfnt,5), ET_comp(pfnt,3), ET_mass(pfnt), ET_np(pfnt))
allocate(PET_balance(pfnt,2))
ET_age = 0.0d0
ET_mass = 0.0d0
ET_comp = 0.0d0
ET_np = 0
PET_balance = 0.0D0

allocate(Out_age(pfnt,5), Out_comp(pfnt,3), Out_mass(pfnt,5), Out_np(pfnt))
Out_age = 0.0d0
Out_mass = 0.0d0
Out_comp = 0.0d0
Out_np = 0

! clear out output particles
npout = 0

!! IO control, each value is a timestep interval, e.g. 1= every timestep, 2=every other, 0 = no writing
read(10,*) ipwrite        ! controls an ASCII, .3D particle file not recommended due to poor performance
read(10,*) ibinpntswrite  !  controls VTK, binary output of particle locations and attributes
read(10,*) etwrite        !  controls ASCII ET output
read(10,*) icwrite        ! controls VTK, binary grid based output where particle masses, concentrations,
                          ! ages are mapped to a grid and written every N timesteps


! allocate and assign timesteps
allocate(Time_Next(pfnt))

do kk = 1, pfnt
        Time_Next(kk) = float(kk+outkk-1)*pfdt
end do

Time_first = float(outkk-1)*pfdt

! read in velocity multiplier
read(10,*) V_mult

! do we read in clm evap trans?
read(10,*) clmtrans
! do we read in clm output file?
read(10,*) clmfile

! read in IC number of particles for flux
read(10,*) iflux_p_res

! read in density h2o
read(10,*) denh2o

! read in diffusivity
!moldiff = (1.15e-9)*3600.d0
read(10,*) moldiff

!saving Efract for a later time
!read(10,*) Efract

! fraction of dx/Vx
read(10,*) dtfrac

!wite out log file
write(11,*)
write(11,*) ' Grid Dimensions'
write(11,'("dx:",e12.5)') dx
write(11,'("dy:",e12.5)') dy
write(11,'("dz:",*(e12.5,", "))') dz(1:nz)
write(11,*)

! end of SLIM input
close(10)

BCAST_int(1)=nx !reduce need nx, ny, and nz
BCAST_int(2)=ny
BCAST_int(3)=nz
BCAST_int(4)=np_ic !not needed
BCAST_int(5)=np
BCAST_int(6)=pft1
BCAST_int(7)=pft2
BCAST_int(8)=tout1
BCAST_int(9)=n_cycle !4-9 may not need
BCAST_int(10)=ipwrite
BCAST_int(11)=ibinpntswrite
BCAST_int(12)=etwrite
BCAST_int(13)=icwrite
BCAST_int(14)=iflux_p_res ! 10-14 may not need
BCAST_int(15)=n_constituents
BCAST_int(16)=pfnt

BCAST_float(1)=dx
BCAST_float(2)=dy
BCAST_float(3)=pfdt
BCAST_float(4)=moldiff
BCAST_float(5)=dtfrac

end if

call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_BCAST(BCAST_int,16,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(BCAST_float,5,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

if(rank /= 0) then

    nx            = BCAST_int(1)
    ny            = BCAST_int(2)
    nz            = BCAST_int(3)
    np_ic         = BCAST_int(4)
    np            = BCAST_int(5)
    pft1          = BCAST_int(6)
    pft2          = BCAST_int(7)
    tout1         = BCAST_int(8)
    n_cycle       = BCAST_int(9)
    ipwrite       = BCAST_int(10)
    ibinpntswrite = BCAST_int(11)
    etwrite       = BCAST_int(12)
    icwrite       = BCAST_int(13)
    iflux_p_res   = BCAST_int(14)
    n_constituents = BCAST_int(15)
    pfnt          = BCAST_int(16)

    dx      = BCAST_float(1) 
    dy      = BCAST_float(2) 
    pfdt    = BCAST_float(3)  
    moldiff = BCAST_float(4) 
    dtfrac  = BCAST_float(5)  
    
    allocate(dz(nz))  
    
end if

call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_BCAST(dz,nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!=============================================================================
! grid +1 variables
nnx=nx+1
nny=ny+1
nnz=nz+1

!allocate arrays
allocate(Vx(nnx,ny,nz), Vy(nx,nny,nz), Vz(nx,ny,nnz))
allocate(Saturation(nx,ny,nz), Porosity(nx,ny,nz),EvapTrans(nx,ny,nz))
allocate(C_p(n_constituents,nx,ny,nz))

!Intialize everything to Zero
Vx = 0.0d0
Vy = 0.0d0
Vz = 0.0d0

Saturation = 0.0D0
Porosity = 0.0D0
EvapTrans = 0.0d0
C_p = 0.0d0


if (rank == 0) then

    write(11,*) 'Timestepping Information'
    write(11,'("ParFlow delta-T, pfdt:",e12.5)') pfdt
    write(11,'("ParFlow timesteps, pfnt:",i12)') pfnt
    write(11,'("ParFlow start step, pft1:",i12)') pft1
    write(11,'("ParFlow end step, pft2:",i12)') pft2
    write(11,'("Output step start:",i12)') outkk
    write(11,'("Time loops, cycles, n_cycle:",i12)') n_cycle
    write(11,'("Total time steps:",i12)') pfnt

    write(11,*)
    write(11,*) 'V mult: ',V_mult,' for forward/backward particle tracking'
    write(11,*) 'CLM Trans: ',clmtrans,' adds / removes particles based on LSM fluxes'
    write(11,*)
    write(11,*) 'Physical Constants'
    write(11,*) 'denh2o: ',denh2o,' M/L^3'
    write(11,*) 'Molecular Diffusivity: ',moldiff,' '
    !write(11,*) 'Fractionation: ',Efract,' '
    write(11,*)
    write(11,*) 'Numerical Stability Information'
    write(11,'("dtfrac: ",e12.5," fraction of dx/Vx")') dtfrac
    
    call system_clock(T2)
    IO_time_read = IO_time_read + (T2-T1)
    
end if

call MPI_Barrier(MPI_COMM_WORLD,ierr)
!=============================================================================
!! set up domain boundaries
Xmin = 0.0d0
Ymin = 0.0d0
Zmin = 0.0d0
Xmax = float(nx)*dx
Ymax = float(ny)*dy
Zmax = 0.0d0
do k = 1, nz
        Zmax = Zmax + dz(k)
end do

call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(rank == 0)then

DEM = 0.0d0
! read in DEM
if (DEMname /= '') then
 fname = trim(adjustl(DEMname))
 call pfb_read(DEM,fname,nx,ny,nztemp)
end if ! DEM

!! hard wire DEM  @RMM, to do, need to make this input
!do i = 1, nx
!  do j = 1, ny
!DEM(i,j) = 0.0D0 + float(i)*dx*0.05
!end do
!end do


write(11,*)
write(11,*) '## Domain Info'
write(11,'("Xmin:",e12.5," Xmax:",e12.5)') Xmin, Xmax
write(11,'("Ymin:",e12.5," Ymax:",e12.5)') Ymin, Ymax
write(11,'("Zmin:",e12.5," Zmax:",e12.5)') Zmin, Zmax

!! DEM set to zero but will be read in as input

!! Set up grid locations for file output
npnts=nnx*nny*nnz
ncell=nx*ny*nz

allocate(Pnts(npnts,3))
Pnts=0
m=1

! Need the maximum height of the model and elevation locations
Z = 0.0d0
Zt(0) = 0.0D0
do ik = 1, nz
Z = Z + dz(ik)
Zt(ik) = Z
!print*, Z, dz(ik), Zt(ik), ik
end do
maxZ=Z

!! candidate loops for OpenMP
do k=1,nnz
 do j=1,nny
  do i=1,nnx
   Pnts(m,1)=DBLE(i-1)*dx
   Pnts(m,2)=DBLE(j-1)*dy
! This is a simple way of handling the maximum edges
   if (i <= nx) then
   ii=i
   else
   ii=nx
   endif
   if (j <= ny) then
   jj=j
   else
   jj=ny
   endif
   ! This step translates the DEM
   ! The specified initial heights in the pfb (z1) are ignored and the
   !  offset is computed based on the model thickness
   Pnts(m,3)=(DEM(ii,jj)-maxZ)+Zt(k-1)
   m=m+1
  end do
 end do
end do


! Read porosity values from ParFlow .pfb file
fname=trim(adjustl(pname))//'.out.porosity.pfb'
call pfb_read(Porosity,fname,nx,ny,nz)

! Read the in initial Saturation from ParFlow
kk = 0
pfkk=pft1-1
!print*, pfkk
write(filenum,'(i5.5)') pfkk
fname=trim(adjustl(pname))//'.out.satur.'//trim(adjustl(filenum))//'.pfb'
call pfb_read(Saturation,fname,nx,ny,nz)

! Intialize random seed
ir = -3333

!! Define initial particles' locations and mass
!!
if (np_ic > 0)  then
np_active = 0

PInLoc=0.0d0
!call srand(333)
do i = 1, nx
do j = 1, ny
do k = 1, nz
  if (np_active < np) then   ! check if we have particles left
  if (Saturation(i,j,k) > 0.0) then ! check if we are in the active domain
  do ij = 1, np_ic
  np_active = np_active + 1
  ii = np_active
  ! assign X, Y, Z locations randomly to each cell
  ! assign X, Y, Z locations randomly to each cell
 P(ii,1) = float(i-1)*dx  +ran1(ir)*dx
 PInLoc(ii,1) = P(ii,1)
 P(ii,2) = float(j-1)*dy  +ran1(ir)*dy
 PInLoc(ii,2) = P(ii,2)
  Z = 0.0d0
  do ik = 1, k
  Z = Z + dz(ik)
  end do

  P(ii,3) = Z -dz(k)*ran1(ir)
  PInLoc(ii,3) = P(ii,3)

        ! assign mass of particle by the volume of the cells
        ! and the water contained in that cell
        P(ii,6) = dx*dy*dz(k)*(Porosity(i,j,k)  &
                 *Saturation(i,j,k))*denh2o*(1.0d0/float(np_ic))
        P(ii,7) = 1.0d0
        P(ii,8) = 1.0d0
        ! set up intial concentrations
        C(1,i,j,k) = C(1,i,j,k) + P(ii,8)*P(ii,6) /  &
        (dx*dy*dz(k)*(Porosity(i,j,k)        &
         *Saturation(i,j,k)))
        C(2,i,j,k) = C(2,i,j,k) + P(ii,8)*P(ii,4)*P(ii,6)
        C(4,i,j,k) = C(4,i,j,k) + P(ii,8)*P(ii,7)*P(ii,6)
        C(3,i,j,k) = C(3,i,j,k) + P(ii,8)*P(ii,6)
end do   ! particles per cell
end if   !  active domain
else
  write(11,*) ' **Warning IC input but no paricles left'
  write(11,*) ' **Exiting code gracefully writing restart '
  goto 9090

end if
end do ! i
end do ! j
end do ! k

!! if np_ic = -1 then we read a restart file
else if (np_ic == -1) then
  write(11,*) 'Reading particle restart File:',trim(runname)//'_particle_restart.bin'
  ! read in full particle array as binary restart file, should name change?,
  ! potential overwrite confusion
  open(116,file=trim(runname)//'_particle_restart.bin', FORM='unformatted',  &
      access='stream')
  read(116) np_active
  if (np_active < np) then   ! check if we have particles left
  !do ii = 1, np_active
            read(116)  P(1:np_active,1:10)
  !end do !11
  close(116)
  else
    write(11,*) ' **Warning restart IC input but no paricles left'
    write(11,*) ' **Exiting code *not* (over)writing restart '
  close(116)
    stop
  end if
end if
flush(11)

!! Define initial particles' locations by surface water
!!
if (np_ic < -1)  then
np_active = 0

PInLoc=0.0d0
!call srand(333)
ir = -3333
do i = 1, nx
do j = 1, ny
k = nz
  if (np_active < np) then   ! check if we have particles left
  !if (saturation(i,j,k) >= 0.95d0)  then
  do ij = 1, abs(np_ic)
  np_active = np_active + 1
  ii = np_active
  ! assign X, Y, Z locations randomly to each cell
  ! assign X, Y, Z locations randomly to each cell
 P(ii,1) = float(i-1)*dx  +ran1(ir)*dx
 PInLoc(ii,1) = P(ii,1)
 P(ii,2) = float(j-1)*dy  +ran1(ir)*dy
 PInLoc(ii,2) = P(ii,2)
  Z = 0.0d0
  do ik = 1, k
  Z = Z + dz(ik)
  end do

  P(ii,3) = Z !-dz(k)*ran1(ir)
  PInLoc(ii,3) = P(ii,3)

        ! assign mass of particle by the volume of the cells
        ! and the water contained in that cell
        P(ii,6) = dx*dy*dz(k)*(Porosity(i,j,k)  &
                 *Saturation(i,j,k))*denh2o*(1.0d0/float(np_ic))
        P(ii,7) = 1.0d0
        P(ii,8) = 1.0d0
        ! set up intial concentrations
        C(1,i,j,k) = C(1,i,j,k) + P(ii,8)*P(ii,6) /  &
        (dx*dy*dz(k)*(Porosity(i,j,k)        &
         *Saturation(i,j,k)))
        C(2,i,j,k) = C(2,i,j,k) + P(ii,8)*P(ii,4)*P(ii,6)
        C(4,i,j,k) = C(4,i,j,k) + P(ii,8)*P(ii,7)*P(ii,6)
        C(3,i,j,k) = C(3,i,j,k) + P(ii,8)*P(ii,6)
end do   ! particles per cell
!end if  !! saturated at top
else
  write(11,*) ' **Warning IC input but no paricles left'
  write(11,*) ' **Exiting code gracefully writing restart '
  goto 9090

end if
end do ! i
end do ! j

end if !! river IC

! Write out intial concentrations
! normalize ages by mass
where (C(3,:,:,:)>0.0)  C(2,:,:,:) = C(2,:,:,:) / C(3,:,:,:)
where (C(3,:,:,:)>0.0)  C(4,:,:,:) = C(4,:,:,:) / C(3,:,:,:)

!! Set up output options for VTK grid output
!icwrite = 1
vtk_file=trim(runname)//'_cgrid'
conc_header(1) = 'Concentration'
conc_header(2) = 'Age'
conc_header(3) = 'Mass'
conc_header(4) = 'Comp'
conc_header(5) = 'Delta'
conc_header(6) = 'ET_Npart'
conc_header(7) = 'ET_Mass'
conc_header(8) = 'ET_Age'
conc_header(9) = 'ET_Comp'

if(icwrite > 0)  &
call vtk_write(Time_first,C,conc_header,nx,ny,nz,0,n_constituents,Pnts,vtk_file)

!! clear out C arrays
C = 0.0D0

!! write timestep header
write(11,*)
write(11,*)
write(11,*) ' **** Transient Simulation Particle Accounting ****'
write(11,*) ' Timestep PFTimestep OutStep    Time     Mean_Age    Mean_Comp   Mean_Mass  Total_Mass    PrecipIn    ETOut  &
              NP_PrecipIn NP_ETOut &
             NP_QOut NP_active_old NP_filtered'
flush(11)

!! open exited partile file and write header
!open(114,file=trim(runname)//'_exited_particles.txt')
open(114,file=trim(runname)//'_exited_particles.bin', FORM='unformatted',  &
    access='stream')

!write(114,*) 'Time X Y Z PTime Mass Comp ExitStatus'
!flush(114)

if(ipwrite < 0) then
! open/create/write the 3D output file which will write particles out each timestemp, very slowly in parallel
open(214,file=trim(runname)//'_total_particle_trace.3D')
write(214,*) 'X Y Z TIME'
end if !! ipwrite < 0

!--------------------------------------------------------------------
! (3) For each timestep, loop over all particles to find and
!     update their new locations
!--------------------------------------------------------------------
! loop over timesteps
pfkk = pft1 - 1
outkk = outkk - 1

end if !rank==0

call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_BCAST(Porosity,nx*ny*nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

do kk = 1, pfnt
!print *, kk
if (rank == 0) then
!! reset ParFlow counter for cycles
if (mod((kk-1),(pft2-pft1+1)) == 0 )  pfkk = pft1 - 1
        call system_clock(T1)
        !adust the file counters
        pfkk = pfkk + 1
        outkk = outkk + 1
        !print*, pfkk

        ! Read the velocities computed by ParFlow
        write(filenum,'(i5.5)') pfkk

        fname=trim(adjustl(pname))//'.out.velx.'//trim(adjustl(filenum))//'.pfb'
        call pfb_read(Vx,fname,nx+1,ny,nz)

        fname=trim(adjustl(pname))//'.out.vely.'//trim(adjustl(filenum))//'.pfb'
        call pfb_read(Vy,fname,nx,ny+1,nz)

        fname=trim(adjustl(pname))//'.out.velz.'//trim(adjustl(filenum))//'.pfb'
        call pfb_read(Vz,fname,nx,ny,nz+1)

        fname=trim(adjustl(pname))//'.out.satur.'//trim(adjustl(filenum))//'.pfb'
        call pfb_read(Saturation,fname,nx,ny,nz)

        if (clmtrans) then
        ! Read in the Evap_Trans
        fname=trim(adjustl(pname))//'.out.evaptrans.'//trim(adjustl(filenum))//'.pfb'
        call pfb_read(EvapTrans,fname,nx,ny,nz)
        ! check if we read full CLM output file
        if (clmfile) then
         !Read in CLM output file @RMM to do make this input
        fname=trim(adjustl(pname))//'.out.clm_output.'//trim(adjustl(filenum))//'.C.pfb'
        call pfb_read(CLMvars,fname,nx,ny,nzclm)
        end if
        end if

        call system_clock(T2)

        IO_time_read = IO_time_read + (T2-T1)

        ! Determine whether to perform forward or backward patricle tracking
        Vx = Vx * V_mult
        Vy = Vy * V_mult
        Vz = Vz * V_mult
        i_added_particles = 0
        ! Add particles if P-ET > 0
        if (clmtrans) then   !check if this is our mode of operation, read in the ParFlow evap trans file
                             ! normally generated by CLM but not exclusively, to assign new particles to any
                             ! additional water fluxes (rain, snow, irrigation water) with the mass of each
                             ! particle assigned to be the mass of the NEW water
        ! check overall if we are out of particles (we do this twice once for speed, again for array)
        ! generally if we are out of particles the simulation isn't valid, but we just warn the user
        if (np_active < np) then
        ! loop over entire domain to check each cell, if the cell flux is a recharge we add particles
        ! and sum the flux, if negative we sum ET flux
        ! this could be parallelized but thread race concerns when summing total fluxes
        ! and a scheme to place particles accurately in parallel makes me think serial is
        ! still more efficient
        do i = 1, nx
        do j = 1, ny
        do k = 1, nz
        if (EvapTrans(i,j,k)> 0.0d0) then
        ! sum water inputs in PET 1 = P, 2 = ET, kk= PF timestep
        ! units of ([T]*[1/T]*[L^3])/[M/L^3] gives Mass of water input
        PET_balance(kk,1) = PET_balance(kk,1) &
                            + pfdt*EvapTrans(i,j,k)*dx*dy*dz(k)*denh2o
        do ji = 1, iflux_p_res
        if (np_active < np) then   ! check if we have particles left
        np_active = np_active + 1
        ii = np_active             ! increase total number of particles
        i_added_particles = i_added_particles + 1   ! increase particle counter for accounting
        ! assign X, Y locations randomly to recharge cell
        P(ii,1) = float(i-1)*dx  +ran1(ir)*dx
        PInLoc(ii,1) = P(ii,1)
        P(ii,2) = float(j-1)*dy  +ran1(ir)*dy
        PInLoc(ii,2) = P(ii,2)
        Z = 0.0d0
        do ik = 1, k
        Z = Z + dz(ik)
        end do
        ! Z location is fixed
        P(ii,3) = Z -dz(k)*0.5d0 !  *ran1(ir)
        PInLoc(ii,3) = P(ii,3)

        ! assign zero time and flux of water
        ! time is assigned randomly over the recharge time to represent flux over the
        ! PF DT unless we are running SS then have all particles start at the same time
        P(ii,4) = 0.0d0
        if (iflux_p_res >= 0) P(ii,4) = 0.0d0 +ran1(ir)*pfdt
        P(ii,5) = 0.0d0
        ! mass of water flux into the cell divided up among the particles assigned to that cell
        !P(ii,6) = (1.0d0/float(iflux_p_res))   &
          !        *P(ii,4)*EvapTrans(i,j,k)*dx*dy*dz(k)*denh2o  !! units of ([T]*[1/T]*[L^3])/[M/L^3] gives Mass
        P(ii,6) = (1.0d0/float(abs(iflux_p_res)))   &
                    *pfdt*EvapTrans(i,j,k)*dx*dy*dz(k)*denh2o  !! units of ([T]*[1/T]*[L^3])/[M/L^3] gives Mass
        !! check if input is rain or snowmelt
        if(CLMvars(i,j,11) > 0.0) then !this is snowmelt
        P(ii,7) = 3.0d0 ! Snow composition
        else
        P(ii,7) = 2.d0 ! Rainfall composition
        end if
        P(ii,8) = 1.0d0   ! make particle active
        P(ii,9) = 1.0d0
        P(ii,10) = 0.0d0  ! Particle hasn't exited domain

        else
        write(11,*) ' **Warning rainfall input but no paricles left'
        write(11,*) ' **Exiting code gracefully writing restart '
        goto 9090
        end if  !! do we have particles left?
        end do !! for flux particle resolution

        else !! ET not P
        ! sum water inputs in PET 1 = P, 2 = ET, kk= PF timestep
        ! units of ([T]*[1/T]*[L^3])/[M/L^3] gives Mass of water input
        PET_balance(kk,2) = PET_balance(kk,2) &
                            + pfdt*EvapTrans(i,j,k)*dx*dy*dz(k)*denh2o
        end if  !! end if for P-ET > 0
        end do
        end do
        end do

        end if  !! second particle check to avoid array loop if we are out of particles
        end if  !! end if for clmtrans logical
        !write(11,*) ' Time Step: ',Time_Next(kk),' NP Active:',np_active

        call system_clock(T1)
        
        np_p = floor(float(np_active)/float(pp))
        do i=1,pp
            if(i <= mod(np_active,pp)) then
            np_p(i) = np_p(i) + 1
            end if
        end do
        
end if ! rank == 0

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(np_active,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(np_p,pp,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)   
! message transfer

    !print *, np_active, np_p, rank
        
    if(rank == 0)then

        do i=1,pp-1
           j=sum(np_p(1:i))+1
           k=sum(np_p(1:i+1))
           call MPI_SEND(P(j:k,:),np_p(i+1)*10,MPI_DOUBLE_PRECISION,i,60,MPI_COMM_WORLD,ierr) 
        end do
    else
        allocate(P(np_p(rank+1),10))
        call MPI_RECV(P,np_p(rank+1)*10,MPI_DOUBLE_PRECISION,0,60,MPI_COMM_WORLD,status,ierr)  
    end if
    
    call MPI_Barrier(MPI_COMM_WORLD,ierr)    
    !broadcast
    call MPI_BCAST(Vx,(nx+1)*ny*nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Vy,nx*(ny+1)*nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Vz,nx*ny*(nz+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Saturation,nx*ny*nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(EvapTrans,nx*ny*nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    !call MPI_BCAST(CLMvars,nx*ny*nzclm,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    
    ET_age_p = 0.d0
    ET_mass_p = 0.d0
    ET_comp_p = 0.d0
    ET_np_p = 0
    
    Out_age_p = 0.d0
    Out_mass_p = 0.d0
    Out_comp_p = 0.d0
    Out_np_p = 0 

    C_p = 0.0D0 
! message transfer end
    if (rank == 0) call system_clock(T3)
	
        do ii = 1, np_p(rank+1)
          delta_time = 0.0d0
        !! skip inactive particles still allocated
        !! set random seed for each particle based on timestep and particle number
        ir = -(932117 + sum(np_p(1:rank)) + ii + 100*kk)
        if(P(ii,8) == 1.0) then
                delta_time = P(ii,4) + pfdt
                do while (P(ii,4) < delta_time)

                        ! Find the "adjacent" cell corresponding to the particle's location
                        Ploc(1) = floor(P(ii,1) / dx)
                        Ploc(2) = floor(P(ii,2) / dy)

                        Z = 0.0d0
                        do k = 1, nz
                                Z = Z + dz(k)
                                if (Z >= P(ii,3)) then
                                        Ploc(3) = k - 1
                                        exit
                                end if
                        end do

                ! check to make sure particles are in central part of the domain and if not
                ! apply some boundary condition to them
                !! check if particles are in domain, need to expand this to include better treatment of BC's
                if ((P(ii,1) < Xmin).or.(P(ii,2)<Ymin).or.(P(ii,3)<Zmin).or.  &
                (P(ii,1)>=Xmax).or.(P(ii,2)>=Ymax).or.(P(ii,3)>=(Zmax-dz(nz)))) then

                 ! if outflow at the top add to the outflow age
                !Z = 0.0d0
                !do k = 1, nz
                !Z = Z + dz(k)
                !end do
                if ( (P(ii,3) >= Zmax-(dz(nz)*0.5d0)).and.   &
                (Saturation(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1)  == 1.0).and.  &
                 (Vz(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) > 0.0) ) then
                itime_loc = kk
                if (itime_loc <= 0) itime_loc = 1
                if (itime_loc >= pfnt) itime_loc = pfnt
                !!$OMP ATOMIC
                Out_age_p = Out_age_p + P(ii,4)*P(ii,6)
                !!$OMP ATOMIC
                Out_mass_p = Out_mass_p  + P(ii,6)

                if (P(ii,7) == 1.0) then
                !!$OMP ATOMIC
                Out_comp_p(1) = Out_comp_p(1) + P(ii,6)
                end if
                if (P(ii,7) == 2.0) then
                !!$OMP ATOMIC
                Out_comp_p(2) = Out_comp_p(2) + P(ii,6)
                end if
                if (P(ii,7) == 3.0) then
                !!$OMP ATOMIC
                Out_comp_p(3) = Out_comp_p(3) + P(ii,6)
                end if
                !Out_comp(itime_loc,1) = Out_comp(itime_loc,1) + P(ii,7)*P(ii,6)

                !!$OMP ATOMIC
                Out_np_p = Out_np_p + 1


!               write(22,220) Time_Next(kk), P(ii,1), P(ii,2), P(ii,3), P(ii,4), P(ii,6), P(ii,7)
!    220         FORMAT(7(e12.5))
!                flush(21)
!                !flag particle as inactive
                P(ii,8) = 0.0d0
                 !flag as exiting via Outflow
                P(ii,10) = 1.0d0
                goto 999

                end if
                ! otherwise we just leave it in the domain to reflect
                end if

                        ! Find each particle's factional cell location
                        Clocx = (P(ii,1) - float(Ploc(1))*dx)  / dx
                        Clocy = (P(ii,2) - float(Ploc(2))*dy)  / dy

                        Z = 0.0d0
                        do k = 1, Ploc(3)
                                Z = Z + dz(k)
                        end do
                        Clocz = (P(ii,3) - Z) / dz(Ploc(3) + 1)

                        ! Calculate local particle velocity using linear interpolation,
                        ! converting darcy flux to average linear velocity

                        Vpx = ((1.0d0-Clocx)*Vx(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) &
                              + Vx(Ploc(1)+2,Ploc(2)+1,Ploc(3)+1)*Clocx)   &
                              /(Porosity(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) &
                              *Saturation(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1))

                        Vpy = ((1.0d0-Clocy)*Vy(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) &
                                + Vy(Ploc(1)+1,Ploc(2)+2,Ploc(3)+1)*Clocy) &
                                /(Porosity(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) &
                                *Saturation(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1))

                        Vpz = ((1.0d0-Clocz)*Vz(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) &
                                  + Vz(Ploc(1)+1,Ploc(2)+1,Ploc(3)+2)*Clocz)  &
                                    /(Porosity(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) &
                                  *Saturation(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1))

                        ! calculate particle dt
                        ! check each direction independently
                        advdt = pfdt
                        if (Vpx /= 0.0d0) advdt(1) = dabs(dtfrac*(dx/Vpx))
                        if (Vpy /= 0.0d0) advdt(2) = dabs(dtfrac*(dy/Vpy))
                        if (Vpz /= 0.0d0) advdt(3) = dtfrac*(dz(Ploc(3)+1)/dabs(Vpz))
                        !if (Vpx > 0.0d0) advdt(1) = dabs(((1.0d0-Clocx)*dx)/Vpx)
                        !if (Vpx < 0.0d0) advdt(1) = dabs((Clocx*dx)/Vpx)
                        !if (Vpy > 0.0d0) advdt(2) = dabs(((1.0d0-Clocy)*dy)/Vpy)
                        !if (Vpy < 0.0d0) advdt(2) = dabs((Clocy*dy)/Vpy)
                        !if (Vpz > 0.0d0) advdt(3) = (((1.0d0-Clocz)*dz(Ploc(3)+1))/dabs(Vpz))
                        !if (Vpz < 0.0d0) advdt(3) = ((Clocz*dz(Ploc(3)+1))/dabs(Vpz))

!                        particledt = min(advdt(1)+1.0E-5,advdt(2)+1.0E-5, advdt(3)+1.0E-5, &
!                                  pfdt*dtfrac  ,delta_time-P(ii,4)+1.0E-5)
                        particledt = min(advdt(1),advdt(2), advdt(3), &
                                  pfdt*dtfrac  ,delta_time-P(ii,4))

                        ! calculate Flux in cell and compare it with the ET flux out of the cell
                        if (EvapTrans(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) < 0.0d0)then

                        ! calculate divergence of Darcy flux in the cell
                        !  in X, Y, Z [L^3 / T]
                        local_flux = (Vx(Ploc(1)+2,Ploc(2)+1,Ploc(3)+1) - Vx(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1)) +  &
                                     (Vy(Ploc(1)+1,Ploc(2)+2,Ploc(3)+1) - Vy(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1)) +  &
                                     (Vz(Ploc(1)+1,Ploc(2)+1,Ploc(3)+2) - Vz(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1))

                        ! calculate ET flux volumetrically and compare to
                        et_flux = abs(EvapTrans(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1))*dx*dy*dz(Ploc(3)+1)

                        ! compare total water removed from cell by ET with total water available in cell to arrive at a particle
                        ! probability of being captured by roots
                        ! water volume in cell
                        water_vol = dx*dy*dz(Ploc(3)+1)*(Porosity(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1)  &
                        *Saturation(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1))
                        !  add that amout of mass to ET BT; check if particle is out of mass
                        itime_loc = kk
                        ! extra checking for array bounds to provide future flexibilty
                        ! if DT for output != PFDT
                        if (itime_loc <= 0) itime_loc = 1
                        if (itime_loc >= pfnt) itime_loc = pfnt
                        Zr = ran1(ir)
                        if (Zr < ((et_flux*particledt)/water_vol)) then   ! check if particle is 'captured' by the roots
                        !  this section made atomic since it could inovlve a data race
                        !  that is, each thread can only update the ET arrays one at a time
                        !!$OMP ATOMIC
                        ET_age_p = ET_age_p + P(ii,4)*P(ii,6)  ! mass weighted age
                        !!$OMP ATOMIC
                        ET_mass_p = ET_mass_p  +  P(ii,6)  ! particle mass added to ET

                        !ET_comp(itime_loc,1) = ET_comp(itime_loc,1) + P(ii,7)*P(ii,6)  !mass weighted contribution
                        if (P(ii,7) == 1.0) then
                        !!$OMP ATOMIC
                        ET_comp_p(1) = ET_comp_p(1) + P(ii,6)
                        end if
                        if (P(ii,7) == 2.0) then
                        !!$OMP ATOMIC
                        ET_comp_p(2) = ET_comp_p(2) + P(ii,6)
                        end if
                        if (P(ii,7) == 3.0) then
                        !!$OMP ATOMIC
                        ET_comp_p(3) = ET_comp_p(3) + P(ii,6)
                        end if

                        !!$OMP ATOMIC
                        ET_np_p = ET_np_p + 1   ! track number of particles

                        !outputting spatially distributed ET information
                        !!$OMP Atomic
                        C_p(6,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) = C_p(6,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) + 1 ! Number of ET particles
                        !!$OMP Atomic
                        C_p(7,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) = C_p(7,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) +  P(ii,6) ! particle mass added to ET
                        !!$OMP Atomic
                        C_p(8,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) = C_p(8,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) + P(ii,4)*P(ii,6)  ! mass weighted age
                        !!$OMP Atomic
                        C_p(9,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) = C_p(9,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) + P(ii,7)*P(ii,6)  ! mass weighted contribution

                        !  now remove particle from domain
                            P(ii,8) = 0.0d0
                            !flag as exiting via ET
                           P(ii,10) = 2.0d0

                            goto 999
                        end if
                        end if

                        ! Advect particle to new location using Euler advection until next time
                        P(ii,1) = P(ii,1) + particledt * Vpx
                        P(ii,2) = P(ii,2) + particledt * Vpy
                        P(ii,3) = P(ii,3) + particledt * Vpz
                        P(ii,4) = P(ii,4) + particledt

                        ! Molecular Diffusion
                        if (moldiff > 0.0d0) then
                        z1 = 2.d0*DSQRT(3.0D0)*(ran1(ir)-0.5D0)
                        z2 = 2.d0*DSQRT(3.0D0)*(ran1(ir)-0.5D0)
                        z3 = 2.d0*DSQRT(3.0D0)*(ran1(ir)-0.5D0)

                        P(ii,1) = P(ii,1) + z1 * DSQRT(moldiff*2.0D0*particledt)
                        P(ii,2) = P(ii,2) + z2 * DSQRT(moldiff*2.0D0*particledt)
                        P(ii,3) = P(ii,3) + z3 * DSQRT(moldiff*2.0D0*particledt)
                        end if

!!  placeholder for other interactions; potentially added later
!!
!!                        if (Ploc(3) == nz-1)  P(ii,9) = P(ii,9) -Efract*particledt*CLMvars(Ploc(1)+1,Ploc(2)+1,7)
!!
                        ! place to track saturated / groundwater time if needed
                        if(Saturation(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) == 1.0) P(ii,5) = P(ii,5) + particledt
                        ! simple reflection boundary
                        if (P(ii,3) >=Zmax) P(ii,3) = Zmax- (P(ii,3) - Zmax)
                        if (P(ii,1) >=Xmax) P(ii,1) = Xmax- (P(ii,1) - Xmax)
                        if (P(ii,2) >=Ymax) P(ii,2) = Ymax- (P(ii,2) - Ymax)
                        if (P(ii,2) <=Ymin) P(ii,2) = Ymin+ (Ymin - P(ii,2))
                        if (P(ii,3) <=Zmin) P(ii,3) = Zmin+ (Zmin - P(ii,3) )
                        if (P(ii,1) <=Xmin) P(ii,1) = Xmin+ (Xmin - P(ii,1) )

                        !! write all active particles at concentration in ASCII VisIT 3D file format continuously
                        !! as noted above, this option is very slow compared to VTK binary output
                        !if (ipwrite < 0) then
                        !  if (P(ii,8) == 1.0d0) write(214,61) P(ii,1), P(ii,2), P(ii,3), P(ii,4)
                        !flush(214)
                        !end if !! ipwrite


                end do  ! end of do-while loop for particle time to next time
        999 continue   ! where we go if the particle is out of bounds

                                !! concentration routine
                                ! Find the "adjacent" "cell corresponding to the particle's location
                                Ploc(1) = floor(P(ii,1) / dx)
                                Ploc(2) = floor(P(ii,2) / dy)
                                Z = 0.0d0
                                do k = 1, nz
                                        Z = Z + dz(k)
                                        if (Z >= P(ii,3)) then
                                                Ploc(3) = k - 1
                                                exit
                                        end if
                                end do
                                !!$OMP Atomic
                        C_p(1,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) = C_p(1,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) + P(ii,8)*P(ii,6) /  &
                        (dx*dy*dz(Ploc(3)+1)*(Porosity(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1)        &
                         *Saturation(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1)))
                        !!$OMP Atomic
                        C_p(2,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) = C_p(2,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) + P(ii,8)*P(ii,4)*P(ii,6)
                        !!$OMP Atomic
                        C_p(4,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) = C_p(4,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) + P(ii,8)*P(ii,7)*P(ii,6)
                        !!$OMP Atomic
                        C_p(3,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) = C_p(3,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) + P(ii,8)*P(ii,6)
                        !!$OMP Atomic
                        C_p(5,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) = C_p(5,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) + P(ii,8)*P(ii,9)*P(ii,6)
        end if   !! end-if; check if particle is active
        end do !  ii,  end particle loop
    
	if (rank == 0) call system_clock(T4)        
    call MPI_Barrier(MPI_COMM_WORLD,ierr)    
        ! add outflow, ET, and C from different processes
    if(rank == 0)then
        do i=1,pp-1
           j=sum(np_p(1:i))+1
           k=sum(np_p(1:i+1))
           call MPI_RECV(P(j:k,:),np_p(i+1)*10,MPI_DOUBLE_PRECISION,i,61,MPI_COMM_WORLD,status,ierr) 
        end do
    else
        call MPI_SEND(P,np_p(rank+1)*10,MPI_DOUBLE_PRECISION,0,61,MPI_COMM_WORLD,ierr)  
        deallocate(P)
    end if
    
    call MPI_Barrier(MPI_COMM_WORLD,ierr)        
    call MPI_Reduce(ET_age_p,    ET_age(kk,1), 1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_Reduce(ET_comp_p(1),ET_comp(kk,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_Reduce(ET_comp_p(2),ET_comp(kk,2),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_Reduce(ET_comp_p(3),ET_comp(kk,3),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_Reduce(ET_mass_p,   ET_mass(kk),  1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_Reduce(ET_np_p,     ET_np(kk),    1,MPI_INTEGER,         MPI_SUM,0,MPI_COMM_WORLD,ierr)
    
    call MPI_Reduce(Out_age_p,    Out_age(kk,1), 1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_Reduce(Out_comp_p(1),Out_comp(kk,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_Reduce(Out_comp_p(2),Out_comp(kk,2),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_Reduce(Out_comp_p(3),Out_comp(kk,3),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_Reduce(Out_mass_p,   Out_mass(kk,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_Reduce(Out_np_p,     Out_np(kk),    1,MPI_INTEGER,         MPI_SUM,0,MPI_COMM_WORLD,ierr)
    
    call MPI_Reduce(C_p,C,n_constituents*nx*ny*nz,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)   
    
if (rank == 0) then
        call system_clock(T2)
        parallel_time = parallel_time + (T2-T1)

    write(*,"(6I20)") kk, np_active, T2-T1,T3-T1,T4-T3,T2-T4
call system_clock(T1)


!! format statements for particle output
61  FORMAT(4(e12.5))
62  FORMAT(4(e12.5))

write(filenumout,'(i5.5)') outkk


! write all active particles at concentration in ASCII VisIT 3D file format
! as noted above, this option is very slow compared to VTK binary output
if(ipwrite > 0) then
if(mod(kk,ipwrite) == 0)  then
! open/create/write the 3D output file
open(14,file=trim(runname)//'_transient_particle.'//trim(adjustl(filenumout))//'.3D')
write(14,*) 'X Y Z TIME'
do ii = 1, np_active
if (P(ii,8) == 1) write(14,61) P(ii,1), P(ii,2), P(ii,3), P(ii,4)
end do
close(14)
end if
end if

! normalize ages by mass
where (C(3,:,:,:)>0.0)  C(2,:,:,:) = C(2,:,:,:) / C(3,:,:,:)
where (C(3,:,:,:)>0.0)  C(4,:,:,:) = C(4,:,:,:) / C(3,:,:,:)
where (C(3,:,:,:)>0.0)  C(5,:,:,:) = C(5,:,:,:) / C(3,:,:,:)

! normalize ET ages by mass
where (C(7,:,:,:)>0.0)  C(8,:,:,:) = C(8,:,:,:) / C(7,:,:,:)
where (C(7,:,:,:)>0.0)  C(9,:,:,:) = C(9,:,:,:) / C(7,:,:,:)

!Write gridded ET outputs to text files
if(etwrite > 0) then
if (mod(kk,etwrite) == 0) then
! open/create/write the 3D output file
open(14,file=trim(runname)//'_ET_summary.'//trim(adjustl(filenumout))//'.txt')
write(14,*) 'X Y Z ET_npart, ET_mass, ET_age, ET_comp, EvapTrans_Rate, Saturation, Porosity'
do i = 1, nx
do j = 1, ny
do k = 1, nz
  if (EvapTrans(i,j,k) < 0.0d0)   &
  write(14,'(3(i6), 7(e13.5))')  i, j, k, C(6,i,j,k), C(7,i,j,k), C(8,i,j,k), &
       C(9,i,j,k), EvapTrans(i,j,k), Saturation(i,j,k), Porosity(i,j,k)
end do
end do
end do
close(14)
end if
end if

! write grid based values ("concentrations")
vtk_file=trim(runname)//'_cgrid'
if(icwrite > 0)  then
if(mod(kk,icwrite) == 0)  &
call vtk_write(Time_Next(kk),C,conc_header,nx,ny,nz,outkk,n_constituents,Pnts,vtk_file)
end if
! write binary particle locations and attributes
vtk_file=trim(runname)//'_pnts'
if(ibinpntswrite > 0)  then
if(mod(kk,ibinpntswrite) == 0)  &
call vtk_write_points(P,np_active,np,outkk,vtk_file, dx, dy,nx,ny, maxZ,dem)
end if
!! reset C
C = 0.0D0
!! reset ET_grid
!!ET_grid = 0.0D0
call system_clock(T2)
IO_time_write = IO_time_write + (T2-T1)

call system_clock(T1)
! open exited particle file
!open(114,file=trim(runname)//'_exited_particle.'//trim(adjustl(filenum))//'.txt')
!write(114,*) 'TIMESTEP X Y Z PTIME COMP MASS ExitSTATUS'
!! sort particles to move inactive ones to the end and active ones up
np_active2 = np_active

! !$OMP PARALLEL SHARED(P, np_active, np_active2) &
! !$OMP& Default(none)
!
mean_age = 0.0d0
mean_comp = 0.0d0
total_mass = 0.0D0
mean_mass = 0.0d0

! ! loop over active particles
! !$OMP DO
do ii = 1, np_active
  !! check if particle is inactive
  if (P(ii,8) == 0.0) then
  ! exchange with the last particle and write out exited particles to file
  !write(114,'(6(e13.5),2(i4))') Time_Next(kk), P(ii,1), P(ii,2), P(ii,3), P(ii,4), P(ii,6), int(P(ii,7)), int(P(ii,10))
  !write(114) sngl(Time_Next(kk)), sngl(P(ii,1)), sngl(P(ii,2)), sngl(P(ii,3)), &
  !        sngl(P(ii,4)), sngl(P(ii,6)), sngl(P(ii,7)), sngl(P(ii,10))
          write(114) Time_Next(kk), P(ii,1), P(ii,2), P(ii,3), &
                  P(ii,4), P(ii,6), P(ii,7), P(ii,10)
          npout = npout + 1

!  !$OMP CRITICAL
  P(ii,:) = P(np_active2,:)
  np_active2 = np_active2 -1
!  !$OMP END CRITICAL
else
! increment mean age, composition and mass
mean_age = mean_age + P(ii,4)*P(ii,6)
mean_comp = mean_comp + P(ii,7)*P(ii,6)
total_mass = total_mass + P(ii,6)
end if
  ! if we have looped all the way through our active particles exit
  if (ii > np_active2) exit
end do ! particles
! !$OMP END DO NOWAIT
! !$OMP END PARALLEL
flush(114)
!close(114)
call system_clock(T2)
sort_time = sort_time + (T2-T1)

! check if we have any active particles
if (total_mass > 0.0d0)  then
  mean_age =  mean_age / total_mass
  mean_comp = mean_comp / total_mass
  mean_mass = total_mass / float(np_active2)
end if

! write out summary of mass, age, particles for this timestep
  write(11,'(3(i10),3(f12.5),4(1x,e12.5,1x),3(i8),2(i12))') kk, pfkk, outkk, Time_Next(kk), mean_age , mean_comp, mean_mass, &
                                          total_mass,  PET_balance(kk,1), PET_balance(kk,2), &
                                          i_added_particles,  &
                                          ET_np(kk), Out_np(kk), np_active,np_active2
  flush(11)

np_active = np_active2

end if

call MPI_Barrier(MPI_COMM_WORLD,ierr)

end do !! timesteps and cycles

! Close 3D file
!close(12)
! close particle ET and outflow files
!close(21)
!close(22)

9090 continue  !! continue statement for running out of particles when assigning precip flux.
               !!  code exits gracefully (writes files and possibly a restart file so the user can
               !!  re-run the simulation)
               
if (rank == 0) then

call system_clock(T1)
! Create/open/write the final particles' locations and residence time, @RMM, updated to make this binary and
! make this contain ALL particle attributes to act as a restart file
!open(13,file=trim(runname)//'_endparticle.txt')
!write(13,*) 'NP X Y Z TIME'
!do ii = 1, np_active
!        write(13,63) ii, P(ii,1), P(ii,2), P(ii,3), P(ii,4), P(ii,5), P(ii,6), PInLoc(ii,1), PInLoc(ii,2), PInLoc(ii,3)
!        63  FORMAT(i10,9(e12.5))
!end do
!flush(13)
!! close end particle file
!close(13)

! Write out full particle array as binary restart file
open(116,file=trim(runname)//'_particle_restart.bin', FORM='unformatted',  &
    access='stream')
write(116) np_active
!do ii = 1, np_active
          write(116)  P(1:np_active,1:10)
!end do !11
close(116)

! close output file
close(114)

! Create/open/write the final particles' locations and residence time
! @RMM, superceded by option for 3D file each timestep, this file
! should be the same as the last 3D file written
!open(13,file=trim(runname)//'_endparticle.3D')
!write(13,*) 'X Y Z TIME'
!do ii = 1, np_active
!        write(13,65)  P(ii,1), P(ii,2), P(ii,3), P(ii,4)
!        65  FORMAT(4(e12.5))
!end do
!flush(13)
!! close end particle file
!close(13)
!! write ET files
!
open(13,file=trim(runname)//'_ET_output.txt')
write(13,*) 'TIME ET_age ET_comp1 ET_comp2 ET_comp3 ET_mass ET_Np'
do ii = 1, pfnt
if (ET_mass(ii) > 0 ) then
ET_age(ii,1) = ET_age(ii,1)/(ET_mass(ii))
ET_comp(ii,1) = ET_comp(ii,1)/(ET_mass(ii))
ET_comp(ii,2) = ET_comp(ii,2)/(ET_mass(ii))
ET_comp(ii,3) = ET_comp(ii,3)/(ET_mass(ii))
end if

!if (ET_np(ii) > 0 ) then
!ET_mass(ii,1) = ET_mass(ii,:)   !/(ET_np(ii))
!end if
!write(13,'(6(e12.5),i12)') float(ii)*ET_dt, ET_age(ii,1), ET_comp(ii,1), &
!                           ET_mass(ii,1), ET_mass(ii,2),ET_mass(ii,3), ET_np(ii)
write(13,'(6(e12.5),i12)') float(ii+tout1-1)*ET_dt, ET_age(ii,1), ET_comp(ii,1), &
                            ET_comp(ii,2), ET_comp(ii,3), ET_mass(ii), ET_np(ii)

64  FORMAT(6(e12.5),i12)
end do
flush(13)
! close ET file
close(13)

!! write Outflow
!
open(13,file=trim(runname)//'_flow_output.txt')
write(13,*) 'TIME Out_age Out_comp1 outcomp2 outcomp3 Out_mass Out_NP'
do ii = 1, pfnt
if (Out_mass(ii,1) > 0 ) then
Out_age(ii,:) = Out_age(ii,:)/(Out_mass(ii,1))
Out_comp(ii,1) = Out_comp(ii,1)/(Out_mass(ii,1))
Out_comp(ii,2) = Out_comp(ii,2)/(Out_mass(ii,1))
Out_comp(ii,3) = Out_comp(ii,3)/(Out_mass(ii,1))
!Out_mass(ii,:) = Out_mass(ii,:)/(Out_mass(ii,:))
end if
!write(13,64) float(ii)*ET_dt, Out_age(ii,1), Out_comp(ii,1), Out_mass(ii,1), Out_np(ii)
write(13,64) float(ii+tout1-1)*ET_dt, Out_age(ii,1), Out_comp(ii,1), &
               Out_comp(ii,2), Out_comp(ii,3), Out_mass(ii,1), Out_np(ii)

if(ipwrite < 0) close(214)

end do
flush(13)
! close ET file
close(13)

!! write P-ET water balance
!
open(13,file=trim(runname)//'_PET_balance.txt')
write(13,*) 'TIME P[kg] ET[kg]'
do ii = 1, pfnt
!write(13,'(3(e12.5,2x))') float(ii)*ET_dt, PET_balance(ii,1), PET_balance(ii,2)
write(13,'(3(e12.5,2x))') float(ii+tout1-1)*ET_dt, PET_balance(ii,1), PET_balance(ii,2)
end do
flush(13)
! close ET file
close(13)
call system_clock(T2)
IO_time_write = IO_time_write + (T2-T1)



        call system_clock(Total_time2)

        Write(11,*)
        Write(11,*) '###  Execution Finished'
        write(11,*)
        write(11,*) npout,' particles exited the domain via outflow or ET.'
        write(11,*)
        Write(11,*) 'Simulation Timing and Profiling:'
        Write(11,'("Total Execution Time (s):",e12.5)') float(Total_time2-Total_time1)/1000.
        Write(11,'("File IO Time Read (s):",e12.5)')float(IO_time_read)/1000.
        Write(11,'("File IO Time Write (s):",e12.5)') float(IO_time_write)/1000.
        Write(11,'("Time Sorting (s):",e12.5)') float(sort_time)/1000.
        Write(11,'("Parallel Particle Time (s):",e12.5)') float(parallel_time)/1000.
        write(11,*)
        ! close the log file
        close(11)
        
end if  !rank == 0

    call MPI_Barrier(MPI_COMM_WORLD,ierr)     
    call MPI_FINALIZE(ierr)
    
   end program EcoSLIM
   !
!-----------------------------------------------------------------------
!     function to generate pseudo random numbers, uniform (0,1)
!-----------------------------------------------------------------------
!
!function rand2(iuu)
!
!real*8 rand2,rssq,randt
!integer m,l,iuu

!
!data m/1048576/
!data l/1027/
!n=l*iuu
!iuu=mod(n,m)
!rand2=float(iuu)/float(m)
!    rand2 = rand(0)
!iiu = iiu + 2
!return
!end
