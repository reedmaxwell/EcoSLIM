!--------------------------------------------------------------------
! SLIM2 is a Lagrangian, particle-tracking model for simulating
! subsurface and surface transport of reactive (such as
! microbial agents and metals) and non-reactive contaminants,
! diagnosing travel times, paths etc., which integrates
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
!--------------------------------------------------------------------
! MAIN FORTRAN CODE
!--------------------------------------------------------------------
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

!--------------------------------------------------------------------
! (1) Define variables
!--------------------------------------------------------------------


program SLIM2

  implicit none

  real*8,allocatable::P(:,:)
  ! P = Particle array [np,attributes]
  ! np = Number of particles
  ! P(np,1) = X [real] coordinate
  ! P(np,2) = Y coordinate
  ! P(np,3) = Z coordinate
  ! P(np,4) = Particle residence time
  ! P(np,5) = Saturated particle residence time
  ! P(np,6) = Particle mass; assigned via rainfall rate (Evap_Trans*density*volume*dT)

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

  real*8,allocatable::Time_Next(:)
  ! Vector of real times at which ParFlow dumps outputs

  real*8,allocatable::dz(:)
  ! Vector of depths in z-direction wrt the bottom of domain
  ! Note: The model can handle variable depths, but assumes
  !       uniform spacing in x and y directions.

  real*8,allocatable::Sx(:,:)  ! Sx: Slopes in x-direction
  real*8,allocatable::Sy(:,:)  ! Sy: Slopes in y-direction



  real*8,allocatable::Saturation(:,:,:)   ! Saturation (read from ParFlow)
  real*8,allocatable::Porosity(:,:,:)     ! Porosity (read from ParFlow)
  real*8,allocatable::EvapTrans(:,:,:)     ! CLM EvapTrans (read from ParFlow, [1/T] units)

  real*8  denh20   ! density of water

  integer nx
  integer ny
  integer nz
  ! Domain number of cells in x, y, and z directions

  integer np_ic, np, np_active
  ! number of particles for intial pulse IC, total, and running active

  integer nt
  ! number of timesteps

  real*8  pfdt
  ! ParFlow delta-T

  integer pfnt
  ! number of ParFlow timesteps


  logical clmtrans
  ! logical for mode of operation wirth CLM, will add particles with P-ET > 0
  ! will remove particles if ET > 0

  real*8 dtfrac
  ! fraction of dx/Vx (as well as dy/Vy and dz/Vz) assuring
  ! numerical stability while advecting a particle to a new
  ! location.

  real*8 Xmin, Xmax, Ymin, Ymax, Zmin, Zmax
  ! Domain boundaries in real coordinates. min values set to zero,
  ! but could be adjusted to match Terrain Following Grid in ParFlow.

  real*8 dx, dy
  ! size of grid cells in x and y directions


  real*8 Xlow, Xhi, Ylow, Yhi, Zlow, Zhi
  ! Particles initial locations i.e., where they are injected
  ! into the domain.

  real*8, allocatable:: ET_age(:,:)
  ! time history of ET, time (1,:) and mass for rain (2,:), snow (3,:)
  real*8  ET_dt
  ! time interval for ET


  !! local variables
  integer kk
  ! Loop counter for the time steps (pfnt)
  integer ii
  ! Loop counter for the number of particles (np)
  integer iii
  ! Parallel loop counter for the number of particles run in threads (ithreads)
  integer ithread, nthreads, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM, iparticle_group
  ! number of parallel threads chosen by user

  integer i,j

  ! Local indices

  character*200 runname, filenum, pname, fname
  ! runname = SLIM runname
  ! filenum = ParFlow file number
  ! pname = ParFlow output runname
  ! fname = Full name of a ParFlow's output

  real*8  x, y
  ! x, y, z dimensions

  real*8 V_mult
  ! Multiplier for forward/backward particle tracking
  ! If V_mult = 1, forward tracking
  ! If V_mult = -1, backward tracking

  real*8 t1, t2, t3, t4, t5
  ! code timing variables

  ! private local variables used only in subroutine, all made
  !  OMP private

  integer Ploc(3)
  ! Particle's location whithin a cell
  integer k
  integer l
  integer ik
  ! Local indices within the code (x, y, z, vector)
  real*8 Clocx, Clocy, Clocz
  ! The fractional location of each particle within it's grid cell
  real*8 Vpx, Vpy, Vpz
  ! Particle velocity in x, y, and z directions
  real*8  z
  ! local coords

  real*8 particledt, delta_time
  ! The time it takes for a particle to displace from
  ! one location to another and the local particle from-to time
  ! for each PF timestep.

  real*8 local_flux, et_flux, water_vol, Zr
  ! The local cell flux convergence
  ! The volumetric ET flux
  ! The availble water volume in a cell
  ! random variable
  ! each are allocated over NP to avoid mixing between threaded particles

  !------------

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
  !       - #14: runname_transient_particle.XXX.3D  (visualizes particles in VisIT, one per timetep)
  !
  ! See README_SLIM2.txt for detailes on the above files.

  ! first timing check
  CALL CPU_TIME(t1)

  ! open SLIM input .txt file
  open (10,file='slimin.txt')

  ! read SLIM run name
  read(10,*) runname

  !print*, runname

  ! read ParFlow run name
  read(10,*) pname
  !print*, pname

  ! open/create/write the output log.txt file. If doesn't exist, it's created.
  open(11,file=trim(runname)//'_log.txt')
  write(11,*) 'SLIM2 Log File'
  write(11,*) 'run name:',trim(runname)

  ! open/create/write the 3D output file
  open(12,file=trim(runname)//'_particle.3D')
  write(12,*) 'X Y Z TIME'

  ! read domain number of cells and number of partciels to be injected
  read(10,*) nx
  read(10,*) ny
  read(10,*) nz
  read(10,*) np_ic
  read(10,*) np

  if (np_ic > np) then
    write(11,*) ' warning NP_IC greater than IC'
    np = np_ic
  end if

  ! write nx, ny, nz, and np in the log file
  write(11,*) 'nx:',nx
  write(11,*) 'ny:',ny
  write(11,*) 'nz:',nz
  write(11,*) 'np IC:',np_ic
  write(11,*) 'np:',np


  ! allocate P, Sx, dz, Vx, Vy, Vz, Saturation, and Porosity arrays
  allocate(P(np,10))   !,Ploc(np,3))  !,Clocx(np), Clocy(np), Clocz(np),Vpx(np), Vpy(np), Vpz(np),particledt(np), delta_time(np),  &
  !! local_flux(np), et_flux(np), water_vol(np), Zr(np))
  allocate(PInLoc(np,10))
  allocate(Sx(nx,ny),Sy(nx,ny))
  allocate(dz(nz))
  allocate(Vx(nx+1,ny,nz), Vy(nx,ny+1,nz), Vz(nx,ny,nz+1))
  allocate(Saturation(nx,ny,nz), Porosity(nx,ny,nz),EvapTrans(nx,ny,nz))

  ! read dx, dy as scalars
  read(10,*) dx
  read(10,*) dy

  ! read dz as an array
  read(10,*) dz(1:nz)

  ! read in (constant for now) ParFlow dt
  read(10,*) pfdt

  ! read in (constant for now) ParFlow nt
  read(10,*) pfnt

  ! allocate and assign timesteps
  allocate(Time_Next(pfnt))

  do kk = 1, pfnt
    Time_Next(kk) = float(kk)*pfdt
  end do

  ! Uncomment the following line if all particles are wanted
  ! to exit the domain -- holds only for steady state case.
  ! For unsteady case, only holds if the maximum of all
  ! particles travel time is less or equal than the ParFlow
  ! running time.

  !Time_Next(pfnt) = Time_Next(pfnt-1) + 1.0E15

  ! read in bounds for the particles initial locations
  read(10,*) Xlow, Xhi
  read(10,*) Ylow, Yhi
  read(10,*) Zlow, Zhi

  ! read in velocity multiplier
  read(10,*) V_mult

  ! read in clm
  read(10,*) clmtrans

  ! fraction of dx/Vx
  read(10,*) dtfrac

  !wite out log file
  write(11,*) 'dx:',dx
  write(11,*) 'dy:',dy
  write(11,*) 'dz:',dz(1:nz)
  write(11,*) 'pfdt:',pfdt
  write(11,*) 'pfnt:',pfnt
  write(11,*) 'Initial Condition Info'
  write(11,*) 'X low:',Xlow,' X high:',Xhi
  write(11,*) 'Y low:',Ylow,' Y high:',Yhi
  write(11,*) 'Z low:',Zlow,' Z high:',Zhi
  write(11,*)
  write(11,*) 'V mult: ',V_mult,' for forward/backward particle tracking'
  write(11,*) 'CLM Trans: ',clmtrans,' adds / removes particles based on LSM fluxes'
  write(11,*) 'dtfrac: ',dtfrac,' fraction of dx/Vx'

  ! end of SLIM input
  close(10)

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

  write(11,*)
  write(11,*) 'Domain Info'
  write(11,*) 'Xmin:',Xmin,' Xmax:',Xmax
  write(11,*) 'Ymin:',Ymin,' Ymax:',Ymax
  write(11,*) 'Zmin:',Zmin,' Zmax:',Zmax

  !! Define initial particles' locations
  P=0.0d0
  PInLoc=0.0d0
  call srand(333)
  do ii = 1, np_ic
    P(ii,1) = Xlow+rand()*(Xhi-Xlow)
    PInLoc(ii,1) = P(ii,1)
    P(ii,2) = Ylow+rand()*(Yhi-Ylow)
    PInLoc(ii,2) = P(ii,2)
    P(ii,3) = Zlow+rand()*(Zhi-Zlow)
    PInLoc(ii,3) = P(ii,3)
  end do
  flush(11)

  np_active = np_ic

  ! Read porosity values from ParFlow .pfb file
  fname=trim(adjustl(pname))//'.out.porosity.pfb'
  !print*, fname
  call pfb_read(Porosity,fname,nx,ny,nz)
  flush(11)

  !--------------------------------------------------------------------
  ! (3) For each timestep, loop over all particles to find and
  !     update their new locations
  !--------------------------------------------------------------------

  ! second timing check
  CALL CPU_TIME(t2)


  ! loop over timesteps
  do kk = 1, pfnt
    CALL CPU_TIME(t4)
    ! Read the velocities computed by ParFlow
    write(filenum,'(i5.5)') kk

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
    end if


    ! Determine whether to perform forward or backward patricle tracking
    Vx = Vx * V_mult
    Vy = Vy * V_mult
    Vz = Vz * V_mult
    CALL CPU_TIME(t5)


    ! Add particles if P-ET > 0
    if (clmtrans) then   !check if this is our mode of operation
      ! check overall if we are out of particles (we do this twice once for speed, again for array)
      if (np_active < np) then
        ! loop over top of domain, need to make this more flexible in input for 2D cases, etc
        do i = 2, nx-1
          do j = 3, 3   !ny-1
            do k = 1, nz
              if (EvapTrans(i,j,k)> 0.0d0)  then
                if (Saturation(i,j,k)< 1.0d0)  then
                  !  check if we have particles left
                  if (np_active < np) then   ! check if we have particles left
                    np_active = np_active + 1
                    ii = np_active
                    ! assign X, Y, Z locations to recharge cell
                    P(ii,1) = float(i)*dx  +rand()*dx
                    PInLoc(ii,1) = P(ii,1)
                    P(ii,2) = float(j)*dy  +rand()*dy
                    PInLoc(ii,2) = P(ii,2)
                    z = 0.0d0
                    do ik = 1, k
                      z = z + dz(ik)
                    end do
                    P(ii,3) = z   !+0.5*dz(k)
                    PInLoc(ii,3) = P(ii,3)
                    ! assign zero time and flux of water
                    P(ii,4) = 0.0d0
                    P(ii,5) = 0.0d0
                    P(ii,6) = (pfdt*EvapTrans(i,j,k)*dx*dy*dz(k))/denh20  !! units of ([T]*[1/T]*[L^3])/[M/L^3] gives Mass
                    !print*, i,j,k,P(ii,1:6),ii,np_active
                  else
                    !write(11,*) ' **Warning rainfall input but no paricles left'
                  end if  !! do we have particles left?
                end if  !! end for Sat < 1
              end if  !! end if for P-ET > 0
            end do
          end do
        end do

      end if  !! second particle check to avoid array loop if we are out of particles
    end if  !! end if for clmtrans logical
    write(11,*) ' Time Step: ',Time_Next(kk),' NP Active:',np_active
    write(11,*) ' File Loading Timing (s): ', t5-t4
    flush(11)

    ! open/create/write the 3D output file
    open(14,file=trim(runname)//'_transient_particle.'//trim(adjustl(filenum))//'.3D')
    write(14,*) 'X Y Z TIME'

    ! loop over active particles
    !!! Begin threaded parallelization
    CALL CPU_TIME(t4)

!!!    !$OMP PARALLEL PRIVATE(ithread,nthreads,iparticle_group)
    !! !     get thread id
    !ithread = OMP_GET_THREAD_NUM()
    !nthreads = OMP_GET_NUM_THREADS()
    !if (ithread ==0) write(11,*) ' Number of threads: ',nthreads

    !iparticle_group = np_active / nthreads

    !write(11,*) ' Thread: ',ithread, ' Stride: ',(ithread)*iparticle_group+1,(ithread+1)*iparticle_group

    !flush(11)

    !!        do ii = (ithread)*iparticle_group+1,(ithread+1)*iparticle_group
    !!!    !$OMP PARALLEL DO

    ! private local variables used only in subroutine
!$OMP PARALLEL DO PRIVATE(Ploc,k,l,ik,Clocx, Clocy, Clocz,Vpx)  &
!$OMP& PRIVATE(Vpy, Vpz,z, particledt, delta_time,local_flux, et_flux) &
!$OMP& PRIVATE(water_vol, Zr,np_active)  

!! !$OMP& SHARED(np, nx, ny, nz, Vx, Vy, Vz, Saturation, Porosity)  &
!! !$OMP& SHARED(dx, dy, dz, EvapTrans, pfdt, denh20, Xmax, Ymax, Zmax) &
!! !$OMP& SHARED(Xmin, Ymin, Zmin, dtfrac) &
!! !$OMP& REDUCTION(+:P)
    do ii = 1, np_active
      delta_time = P(ii,4) + pfdt
      do while (P(ii,4) < delta_time)

        ! Find the "adjacent" "cell corresponding to the particle's location
        Ploc(1) = floor(P(ii,1) / dx)
        Ploc(2) = floor(P(ii,2) / dy)

        z = 0.0d0
        do k = 1, nz
          z = z + dz(k)
          if (z >= P(ii,3)) then
            Ploc(3) = k - 1
            exit
          end if
        end do

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

        Vpx = ((1.0d0-Clocx)*Vx(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) + Vx(Ploc(1)+2,Ploc(2)+1,Ploc(3)+1)*Clocx)  &
        /(Porosity(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1)*Saturation(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1))

        Vpy =  ((1.0d0-Clocy)*Vy(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) + Vy(Ploc(1)+1,Ploc(2)+2,Ploc(3)+1)*Clocy) &
        /(Porosity(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1)*Saturation(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1))

        Vpz = ((1.0d0-Clocz)*Vz(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) + Vz(Ploc(1)+1,Ploc(2)+1,Ploc(3)+2)*Clocz)  &
        /(Porosity(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1)*Saturation(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1))

        ! calculate particle dt
        particledt = min(dabs(dtfrac*(dx/Vpx)),dabs(dtfrac*(dy/Vpy)),dtfrac*(dz(Ploc(3)+1)/dabs(Vpz)),pfdt)

        ! calculate Flux in cell and compare it with the ET flux out of the cell
        if (EvapTrans(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) < 0.0d0)  then

          ! calculate divergence of Darcy flux in the cell
          !  in X, Y, Z [L^3 / T]
          local_flux = (Vx(Ploc(1)+2,Ploc(2)+1,Ploc(3)+1) - Vx(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1)) +  &
          (Vy(Ploc(1)+1,Ploc(2)+2,Ploc(3)+1) - Vy(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1)) +  &
          (Vz(Ploc(1)+1,Ploc(2)+1,Ploc(3)+2) - Vz(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1))

          ! calculate ET flux volumetrically and compare to
          et_flux = EvapTrans(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1)*dx*dy*dz(Ploc(3)+1)
          ! compare total water removed from cell by ET with total water available in cell to arrive at a particle
          ! probability of being captured by roots
          ! water volume in cell
          water_vol = dx*dy*dz(Ploc(3)+1)*   &
          (Porosity(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1)*Saturation(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1))
          Zr = rand()
          if (Zr < ((et_flux*particledt)/water_vol)) then   ! check if particle is 'captured' by the roots
            ! subtrack flux from particle
            P(ii,6) = P(ii,6) - water_vol*particledt*denh20
            !  add that amout of mass to ET BT; check if particle is out of mass

          end if
        end if

        ! Advect particle to new location using Euler advection until next time

        ! Update particle location
        P(ii,1) = P(ii,1) + particledt * Vpx
        P(ii,2) = P(ii,2) + particledt * Vpy
        P(ii,3) = P(ii,3) + particledt * Vpz
        P(ii,4) = P(ii,4) + particledt

        !
        if(Saturation(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) == 1.0) P(ii,5) = P(ii,5) + particledt

        ! Check if the particle is still in domain. If not, go to the next particle
        if ((P(ii,1) < Xmin).or.(P(ii,2)<Ymin).or.(P(ii,3)<Zmin).or. &
        (P(ii,1)>Xmax).or.(P(ii,2)>Ymax).or.(P(ii,3)>Zmax)) then
        P(ii,6) = 1
        ! Write out the particle's new location X, Y, Z and its time
        ! in the 3D file
        write(12,61) P(ii,1), P(ii,2), P(ii,3), P(ii,4)
        61  FORMAT(4(e12.5))
        flush(12)
        exit
      else
        P(ii,6) = 0
        ! Write out the particle's new location X, Y, Z and its time
        ! in the 3D file
        write(12,62)  P(ii,1), P(ii,2), P(ii,3), P(ii,4)
        62  FORMAT(4(e12.5))
        flush(12)
      end if

      write(14,61)  P(ii,1), P(ii,2), P(ii,3), P(ii,4)


    end do  ! end of do-while loop
    end do
!$OMP END PARALLEL DO


  CALL CPU_TIME(t5)
  write(11,*) ' Particle Loop Timing (s): ', t5-t4
  flush(11)

  close(14)

end do !! timesteps

! Close 3D file
close(12)


! Create/open/write the final particles' locations and residence time
open(13,file=trim(runname)//'_endparticle.txt')
write(13,*) 'X Y Z TIME'
do ii = 1, np_active
  write(13,63) ii, P(ii,1), P(ii,2), P(ii,3), P(ii,4), P(ii,5), P(ii,6), PInLoc(ii,1), PInLoc(ii,2), PInLoc(ii,3)
  63  FORMAT(i10,9(e12.5))
end do
flush(13)
! close end particle file
close(13)

! third timing check
CALL CPU_TIME(t3)

Write(11,*)
write(11,*) '**** End of Run ***'
write(11,*) ' Setup time (s):' , t2-t1
write(11,*) ' Particle tims (s):', t3-t2
write(11,*) ' Total run time (s):',  t3-t1

! close the log file
close(11)
end program SLIM2
