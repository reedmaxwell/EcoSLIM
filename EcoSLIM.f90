!--------------------------------------------------------------------
! EcoSLIM is a Lagrangian, particle-tracking model for simulating
! subsurface and surface transport of reactive (such as
! microbial agents and metals) and non-reactive contaminants,
! diagnosing travel times, paths etc., which integrates
! seamlessly with ParFlow.
!
! Developed by: Reed Maxwell-August 2016 (rmaxwell@mines.edu)
!
! Contributors: Mohammad Danesh-Yazdi (danesh@mines.edu)
!               Lindsay Bearup (lbearup@usbr.gov)
!
! Explore https://inside.mines.edu/~rmaxwell/maxwell_software.shtml
! for examples of ParFlow and SLIM (older version of SLIM2) working
! together.
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

program SLIM2

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
        ! P(np,9) = isotopic concentration

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

integer Ploc(3)   
        ! Particle's location whithin a cell

integer nx   
integer ny
integer nz
        ! Domain's number of cells in x, y, and z directions

integer np_ic, np, np_active
        ! number of particles for intial pulse IC, total, and running active

integer nt  
        ! number of timesteps ParFlow 

real*8  pfdt, advdt(3)
        ! ParFlow timesteps, advection DT for each direction used to chose optimal particle tiemstep

integer pfnt 
        ! number of ParFlow timesteps

integer kk  
        ! Loop counter for the time steps (pfnt)

integer ii  
        ! Loop counter for the number of particles (np)
integer iflux_p_res
        ! Number of parricles per cell for flux input
integer i
integer j
integer k
integer l
integer ik
integer itime_loc
        ! Local indices within the code (x, y, z, vector)
integer ir

character*200 runname, filenum, pname, fname
        ! runname = SLIM runname
        ! filenum = ParFlow file number
        ! pname = ParFlow output runname
        ! fname = Full name of a ParFlow's output

real*8 Clocx, Clocy, Clocz    
        ! The fractional location of each particle within it's grid cell

real*8 V_mult        
        ! Multiplier for forward/backward particle tracking
        ! If V_mult = 1, forward tracking
        ! If V_mult = -1, backward tracking

logical clmtrans
        ! logical for mode of operation wirth CLM, will add particles with P-ET > 0
        ! will remove particles if ET > 0

real*8 dtfrac, ran1
        ! fraction of dx/Vx (as well as dy/Vy and dz/Vz) assuring 
        ! numerical stability while advecting a particle to a new
        ! location.

real*8 Xmin, Xmax, Ymin, Ymax, Zmin, Zmax     
        ! Domain boundaries in real coordinates. min values set to zero, 
        ! but could be adjusted to match Terrain Following Grid in ParFlow.

real*8 dx, dy
        ! Domain's number of cells in x and y directions

real*8 Vpx, Vpy, Vpz  
        ! Particle velocity in x, y, and z directions

real*8 particledt, delta_time
        ! The time it takes for a particle to displace from 
        ! one location to another and the local particle from-to time
        ! for each PF timestep

real*8 local_flux, et_flux, water_vol, Zr, z1, z2, z3
        ! The local cell flux convergence
        ! The volumetric ET flux
        ! The availble water volume in a cell
        ! random variable

real*8 Xlow, Xhi, Ylow, Yhi, Zlo, Zhi
        ! Particles initial locations i.e., where they are injected
        ! into the domain.

! density of water (M/L3), molecular diffusion (L2.T), fractionation
real*8 denh2o, moldiff, Efract

real*8, allocatable::ET_age(:,:), ET_mass(:,:), ET_comp(:,:)
real*8, allocatable::Out_age(:,:), Out_mass(:,:), Out_comp(:,:)
integer, allocatable:: ET_np(:), Out_np(:)
! time history of ET, time (1,:) and mass for rain (2,:), snow (3,:)
real*8  ET_dt
        ! time interval for ET

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
write(11,*) 'EcoSLIM Log File'
write(11,*) 'run name:',trim(runname)

! open/create/write the 3D output file
open(12,file=trim(runname)//'_particle.3D')
write(12,*) 'X Y Z TIME'

! open/create/write ET particle output file
open(21,file=trim(runname)//'_ET_particle.txt')
write(21,*) 'TIME X Y Z Age Flux Source'

! open/create/write Outflow particle output file
open(22,file=trim(runname)//'_Outflow_particle.txt')
write(22,*) 'TIME X Y Z Age Flux Source'

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
allocate(P(np,10))
P(1:np,1:6) = 0    ! clear out all particle attributes
P(1:np,7:9) = 1.0  ! make all particles active to start with and original from 1 = GW/IC

allocate(PInLoc(np,3))
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

! set ET DT to ParFlow one and allocate ET arrays accordingly
ET_dt = pfdt
allocate(ET_age(pfnt,5), ET_comp(pfnt,5), ET_mass(pfnt,5), ET_np(pfnt))
ET_age = 0.0d0
ET_mass = 0.0d0
ET_comp = 0.0d0
ET_np = 0

allocate(Out_age(pfnt,5), Out_comp(pfnt,5), Out_mass(pfnt,5), Out_np(pfnt))
Out_age = 0.0d0
Out_mass = 0.0d0
Out_comp = 0.0d0
Out_np = 0

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

! read in clm flux
read(10,*) clmtrans

! read in IC number of particles for flux
read(10,*) iflux_p_res

! read in density h2o
read(10,*) denh2o

!! right now, hard wire moldiff as effective rate from Barnes/Allison 88
moldiff = (1.15e-9)*3600.d0

!! right now, hard wire evap fractionation as effective rate from Barnes/Allison 88
Efract = (1.15e-9)*3600.d0

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
write(11,*) 'denh2o: ',denh2o,' water density'

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

PInLoc=0.0d0
call srand(333)
ir = 3333
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

! loop over timesteps
do kk = 1, pfnt  
       
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

        ! Add particles if P-ET > 0
        if (clmtrans) then   !check if this is our mode of operation
        ! check overall if we are out of particles (we do this twice once for speed, again for array)
        if (np_active < np) then
        ! loop over top of domain, need to make this more flexible in input for 2D cases, etc
        do i = 2, nx-1
        do j = 2, ny-1
        do k = 1, nz
        if (EvapTrans(i,j,k)> 0.0d0)  then
        if (Saturation(i,j,k)< 1.0d0)  then
        do ji = 1, iflux_p_res
        !  check if we have particles left
        if (np_active < np) then   ! check if we have particles left
        np_active = np_active + 1
        ii = np_active
        ! assign X, Y, Z locations to recharge cell
        P(ii,1) = float(i-1)*dx  +rand()*dx
        PInLoc(ii,1) = P(ii,1)
        P(ii,2) = float(j-1)*dy  +rand()*dy
        PInLoc(ii,2) = P(ii,2)
        z = 0.0d0
        do ik = 1, k
        z = z + dz(ik)
        end do
        P(ii,3) = z  -dz(k)*rand()
        PInLoc(ii,3) = P(ii,3)
        !print*, P(ii,1), P(ii,2), P(ii,3), z, k
        ! assign zero time and flux of water
        P(ii,4) = 0.0d0
        P(ii,5) = 0.0d0
        P(ii,6) = (1.0d0/float(iflux_p_res))*(pfdt*EvapTrans(i,j,k)*dx*dy*dz(k))*denh2o  !! units of ([T]*[1/T]*[L^3])/[M/L^3] gives Mass
        P(ii,7) = 2.d0 ! Rainfall composition
        !print*, i,j,k,P(ii,1:6),ii,np_active
        else
        write(11,*) ' **Warning rainfall input but no paricles left'
        end if  !! do we have particles left?
        end do !! for flux particle resolution
        end if  !! end for Sat < 1
        end if  !! end if for P-ET > 0
        end do
        end do
        end do

        end if  !! second particle check to avoid array loop if we are out of particles
        end if  !! end if for clmtrans logical
        write(11,*) ' Time Step: ',Time_Next(kk),' NP Active:',np_active
        flush(11)

    ! open/create/write the 3D output file
    open(14,file=trim(runname)//'_transient_particle.'//trim(adjustl(filenum))//'.3D')
    write(14,*) 'X Y Z TIME'
    flush(14)
        ! loop over active particles
        do ii = 1, np_active
        !! skip inactive particles still allocated
        if(P(ii,8) == 1.0) then
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

                ! check to make sure particels are in central part of the domain and if not
                ! apply some boundary condition to them
                !! check if particles are in domain, need to expand this to include better treatment of BC's
                if ((P(ii,1) < Xmin).or.(P(ii,2)<Ymin).or.(P(ii,3)<Zmin).or.    &
                (P(ii,1)>=Xmax).or.(P(ii,2)>=Ymax).or.(P(ii,3)>=Zmax)) then

                 ! if outflow at the top add to the outflow age
                z = 0.0d0
                do k = 1, nz
                z = z + dz(k)
                end do
                if ( (P(ii,3) >= z-dz(nz)/2.0d0).and.   &
                (Saturation(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1)  == 1.0) ) then
                itime_loc = kk
                if (itime_loc <= 0) itime_loc = 1
                if (itime_loc >= pfnt) itime_loc = pfnt
                Out_age(itime_loc,1) = Out_age(itime_loc,1) + P(ii,4)
                Out_mass(itime_loc,1) = Out_mass(itime_loc,1)  + P(ii,6)
                Out_comp(itime_loc,1) = Out_comp(itime_loc,1) + P(ii,7)
                Out_np(itime_loc) = Out_np(itime_loc) + 1

               write(22,220) Time_Next(kk), P(ii,1), P(ii,2), P(ii,3), P(ii,4), P(ii,6), P(ii,7)
    220         FORMAT(7(e12.5))
                flush(21)
!                !flag particle as inactive
                P(ii,8) = 0
                goto 999

                end if
                end if


!! broken logic here for history; should remove when BC approach is finalized
! Check if the particle is still in domain. If not, go to the next particle
!  and tag this particle as inactive
!if ((P(ii,1) < Xmim).or.(P(ii,2)<Ymin).or.(P(ii,3)<Zmin).or. &
!(P(ii,1)>Xmax).or.(P(ii,2)>Ymax).or.(P(ii,3)>Zmax)) then
!Write out the particle's new location X, Y, Z and its time
!in the 3D file
!write(12,61) P(ii,1), P(ii,2), P(ii,3), P(ii,4)
!61  FORMAT(4(e12.5))
!flush(12)
!exit
!else
!P(ii,8) = 1.0
! Write out the particle's new location X, Y, Z and its time
! in the 3D file
!write(12,62) P(ii,1), P(ii,2), P(ii,3), P(ii,4)
!62  FORMAT(4(e12.5))
!flush(12)
!end if
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
                        ! check each direction independently
                        advdt = pfdt
                        if (Vpx /= 0.0d0) advdt(1) = dabs(dtfrac*(dx/Vpx))
                        if (Vpy /= 0.0d0) advdt(2) = dabs(dtfrac*(dx/Vpy))
                        if (Vpz /= 0.0d0) advdt(3) = dtfrac*(dz(Ploc(3)+1)/dabs(Vpz))
                        particledt = min(advdt(1),advdt(2), advdt(3),pfdt)

!print*, "dt",particledt, pfdt, advdt(1:3)
!print*, "loc",Clocx, Clocy, Clocz
!print*, "Vel",Vx(Ploc(1)+2,Ploc(2)+1,Ploc(3)+1), Vy(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1), Vz(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1)
!print*, "P",P(ii, 1:4)
                        ! calculate Flux in cell and compare it with the ET flux out of the cell
                        if (EvapTrans(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) < 0.0d0)  then
!print*, EvapTrans(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1), Ploc(1)+1,Ploc(2)+1,Ploc(3)+1

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
                        Zr = ran1(ir)
!                        print*, P(ii,6),P(ii,7),et_flux, water_vol, et_flux*particledt*denh2o !Zr,(et_flux*particledt)/water_vol,Ploc(1)+1,Ploc(2)+1,Ploc(3)+1
  !                      if (Zr < ((et_flux*particledt)/water_vol)) then   ! check if particle is 'captured' by the roots

                        !  add that amout of mass to ET BT; check if particle is out of mass
                        itime_loc = kk
!                        print*, itime_loc, P(ii,4), ET_dt
                        if (itime_loc <= 0) itime_loc = 1
                        if (itime_loc >= pfnt) itime_loc = pfnt
                        ET_age(itime_loc,1) = ET_age(itime_loc,1) + P(ii,4)
                        ET_mass(itime_loc,1) = ET_mass(itime_loc,1)  + et_flux*particledt*denh2o  ! P(ii,6)
                        ET_comp(itime_loc,1) = ET_comp(itime_loc,1) + P(ii,7)
                        ET_np(itime_loc) = ET_np(itime_loc) + 1
                        ! subtract flux from particle, remove from domain
                        P(ii,6) = P(ii,6) - et_flux*particledt*denh2o

                        if (P(ii,6) <= 0.0d0) then
                        write(21,220) Time_Next(kk), P(ii,1), P(ii,2), P(ii,3), P(ii,4), et_flux*particledt*denh2o, P(ii,7)
                        flush(21)
                            P(ii,8) = 0.0d0
                            goto 999
                            end if
                        !end if
                        end if

                        ! Advect particle to new location using Euler advection until next time

                        ! Update particle location
                        P(ii,1) = P(ii,1) + particledt * Vpx
                        P(ii,2) = P(ii,2) + particledt * Vpy
                        P(ii,3) = P(ii,3) + particledt * Vpz
                        P(ii,4) = P(ii,4) + particledt

                        ! Molecular Diffusion
                        z1 = 2.d0*DSQRT(3.0D0)*(rand()-0.5D0)
                        z2 = 2.d0*DSQRT(3.0D0)*(rand()-0.5D0)
                        z3 = 2.d0*DSQRT(3.0D0)*(rand()-0.5D0)

                        P(ii,1) = P(ii,1) + z1 * DSQRT(moldiff*2.0D0*particledt)
                        P(ii,2) = P(ii,2) + z2 * DSQRT(moldiff*2.0D0*particledt)
                        P(ii,3) = P(ii,3) + z3 * DSQRT(moldiff*2.0D0*particledt)

!!  Apply fractionation if we are in the top cell
!!
                        if (Ploc(3) == nz-1)  P(ii,9) = P(ii,9) -Efract*particledt*100.
                        ! changes made in Ploc           
                       if(Saturation(Ploc(1)+1,Ploc(2)+1,Ploc(3)+1) == 1.0) P(ii,5) = P(ii,5) + particledt
                        ! simple reflection
                        if (P(ii,3) >=Zmax) P(ii,3) = Zmax- (P(ii,3) - Zmax)
                        if (P(ii,1) >=Xmax) P(ii,1) = Xmax- (P(ii,1) - Xmax)
                        if (P(ii,2) >=Ymax) P(ii,2) = Ymax- (P(ii,2) - Ymax)
                        if (P(ii,2) <=Ymin) P(ii,2) = Ymin+ (Ymin - P(ii,2))
                        if (P(ii,3) <=Zmin) P(ii,3) = Zmin+ (Zmin - P(ii,3) )
                        if (P(ii,1) <=Xmin) P(ii,1) = Xmin+ (Xmin - P(ii,1) )

                    write(14,61) P(ii,1), P(ii,2), P(ii,3), P(ii,9)

                end do  ! end of do-while loop for particle time to next time
        999 continue   ! where we go if the particle is out of bounds
        end if   !! check if particle is active
        ! format statements lying around, should redo the way this is done
        61  FORMAT(4(e12.5))
        62  FORMAT(4(e12.5))
        end do !  ii,  end particle loop

close(14)



end do !! timesteps

! Close 3D file
close(12)
! close particle ET and outflow files
close(21)
close(22)

! Create/open/write the final particles' locations and residence time
open(13,file=trim(runname)//'_endparticle.txt')
write(13,*) 'NP X Y Z TIME'
do ii = 1, np_active
        write(13,63) ii, P(ii,1), P(ii,2), P(ii,3), P(ii,4), P(ii,5), P(ii,6), PInLoc(ii,1), PInLoc(ii,2), PInLoc(ii,3)
        63  FORMAT(i10,9(e12.5))
end do
flush(13)
! close end particle file
close(13)

!! write ET files
!
open(13,file=trim(runname)//'_ET_output.txt')
write(13,*) 'TIME ET_age ET_comp ET_mass  ET_Np'
do ii = 1, pfnt
if (ET_np(ii) > 0 ) then
ET_age(ii,:) = ET_age(ii,:)/float(ET_np(ii))
ET_comp(ii,:) = ET_comp(ii,:)/float(ET_np(ii))
ET_mass(ii,:) = ET_mass(ii,:)/float(ET_np(ii))
end if
write(13,64) float(ii)*ET_dt, ET_age(ii,1), ET_comp(ii,1), ET_mass(ii,1), ET_np(ii)
64  FORMAT(4(e12.5),i12)
end do
flush(13)
! close ET file
close(13)

!! write Outflow
!
open(13,file=trim(runname)//'_flow_output.txt')
write(13,*) 'TIME Out_age Out_comp Out_mass Out_NP'
do ii = 1, pfnt
if (Out_np(ii) > 0 ) then
Out_age(ii,:) = Out_age(ii,:)/float(Out_np(ii))
Out_comp(ii,:) = Out_comp(ii,:)/float(Out_np(ii))
Out_mass(ii,:) = Out_mass(ii,:)/float(Out_np(ii))
end if
write(13,64) float(ii)*ET_dt, Out_age(ii,1), Out_comp(ii,1), Out_mass(ii,1), Out_np(ii)

end do
flush(13)
! close ET file
close(13)
! close the log file
close(11)


end program SLIM2

!
!-----------------------------------------------------------------------
!     function to generate pseudo random numbers, uniform (0,1)
!-----------------------------------------------------------------------
!
function rand2(iuu)
!
real*8 rand2,rssq,randt
integer m,l,iuu

!
data m/1048576/
data l/1027/
n=l*iuu
iuu=mod(n,m)
rand2=float(iuu)/float(m)
!    rand2 = rand(0)
iiu = iiu + 2
return
end

!
FUNCTION ran1(idum)
INTEGER*4 idum,IA,IM,IQ,IR,NTAB,NDIV
REAL*8 ran1,AM,EPS,RNMX
PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,  &
NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
INTEGER j,k,iv(NTAB),iy
SAVE iv,iy
DATA iv /NTAB*0/, iy /0/
if (idum.le.0.or.iy.eq.0) then
idum=max(-idum,1)
do 11 j=NTAB+8,1,-1
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM
if (j.le.NTAB) iv(j)=idum
11      continue
iy=iv(1)
endif
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM

j=1+iy/NDIV
iy=iv(j)
iv(j)=idum
ran1=min(AM*iy,RNMX)
return
END
