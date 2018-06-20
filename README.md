EcoSLIM
=======

**EcoSLIM** private repo for development.  Will be moved to public repo under `/ParFlow` when published

**EcoSLIM** is a Lagrangian, particle-tracking code that simulates advective and diffusive movement of water parcels.  This code can be used to simulate age, diagnose travel times, source water composition and flowpaths.  It integrates seamlessly with **ParFlow-CLM**.

### Development Team
+ Reed Maxwell <rmaxwell@mines.edu>
+ Mohammad Danesh-Yazdi <danesh@mines.edu>
+ Laura Condon <lecondon@syr.edu>
+ Lindsay Bearup <lbearup@usbr.gov>

### Citatation
For more details on the model and if you use EcoSLIM in published work please cite the following reference:  
   *Maxwell, R.M., L.E. Condon, M. Danesh-Yazdi and L.A. Bearup, Exploring source water mixing and transient residence time distributions of outflow and  evapotranspiration with an integrated hydrologic model and Lagrangian particle tracking approach, in review, 2018.*

Building and Running
--------------------
To build **EcoSLIM** simply type `make` in the directory with the main directory where `EcoSLIM.f90` sits
To set the number of parallel threads use either
`export OMP_NUM_THREADS=16` for bash or
`setenv OMP_NUM_THREADS 16` for t/c-shell.

To run you will need to have a completed ParFlow simulation and an
EcoSLIM input file that must be named `slimin.txt` and follow the
format described below. Note that the slim input file does not need to be co-located
with the ParFlow simulation.  

To run simply execute `EcoSLIM.exe` from the directory that contains the
`slimin.txt` input file.

Refer to the **Examples** directory described below for example workflows
running with ParFlow and **EcoSLIM**

`slimin.txt`  Main input file. Includes domain geometry, **ParFlow** timing and input, total number of particles,   initial conditions and information about **CLM**.

Inputs
--------------------
**Mandatory Inputs**
1. `slimin.txt`: EcoSLIM input file (See following section for description)
2. *ParFlow* outputs for at least one time step written as pfbs:
   * Velocity - In x, y, and z directions. Used for particle advection
   * Saturation - Used to determine the volume of water in a cell for determining particle mass
   * Porosity - Used to determine the volume of water in a cell for determining particle mass

**Optional Inputs**
1. DEM for ParFlow simulation: If no DEM is provided all elevations are set to zero. **##FILL in **
2. Evapotranspiration: This is normally written by *CLM* but can also be written by *ParFlow* (*out.evaptrans.filenumber.pfb* if clmtrans=T-- EvapTrans)**##FILL in**
3. CLM single file output (*.out.clm_output.t.C.pfb* if clmfile=T --CLMvars) used to specify rain and slow **##FILL in**

EcoSLIM Parameters
--------------------
The `slimin.txt` file contains all of the settings for an EcoSLIM Simulation.
An example of this file is provided below and other examples are also provided
in the **Examples** folder. Inputs *must* be provided in the oder they appear in the
template shown here.

Here we describe each parameter and what they do:
1. **EcoSLIM Runname (runname):** The runname for all of the file outputs
* **ParFlow Runname (pname):** The runname used for the ParFlow simulations. If the *ParFlow*
outputs are not located in the same directory as the *EcoSLIM* run this should also
include the directory path to the *ParFlow* simulation as shown in the example below.
* **DEM File Name (DEMname):** The name of the DEM file for the *ParFlow* simulations. If this line is left blank all elevations will be set to zero **##Fill in*
* **ParFlow nx (nx):** Number of grid cells in the x-direction for ParFlow domain
* **ParFlow ny (ny):**  Number of grid cells in the y-direction for ParFlow domain
* **ParFlow nz (nz):**  Number of grid cells in the z-direction for ParFlow domain
* **Number of Initial Particles (np_ic):** If a positive integer is provided this will be the number of
particles placed in every grid cell at the start of the simulation. To start from a previous
*EcoSLIM* output set this value to -1. In this case particles will be initialized from the
`runname_particle_restart.bin` file (refer to the outputs section for details on this file).  **##Check**
* **Total Particles (np)**: The total number of particles allowed in the simulation. If the particle count
exceeds this at any point (i.e. through particle addition with initial conditions or rainfall events)
the simulation will exit.
* **ParFlow dx (dx):** *ParFlow* grid cell size in the x-direction
* **ParFlow dy (dy):** *ParFlow* grid cell size in the y-direction
* **ParFlow dz (dz):** *ParFlow* grid cell size in the x-direction. This should be a list separated
by comas that is nz long (refer to example below).
* **ParFlow time step (pfdt):** The time step used for the ParFlow simulation. The time units of this are determined by the *ParFlow* simulation and all *EcoSLIM* outputs will have the same units. Currently EcoSLIM assumes a constant time step so this is not compatible with growth times steps in *ParFlow*
* **Starting ParFlow File Number (pft1):** The file number for the *ParFlow* output to start the *EcoSLIM* simulation from. Note that the initial conditions will be set from the file number before this **##Clarify IC use ##
* **Ending ParFlow File Number (pft2):** The file number for the *ParFlow* output to stop the *EcoSLIM* simulation at.
* **EcoSLIM Output Start Counter (tout1):** This initializes the file numbering for the *EcoSLIM* outputs. If this is set to zero then the first *EcoSLIM* output number will be set to match  starting *ParFlow* file number specified above.
* **Time Sequence Repeat (n_cycle):** If the time sequence repeat is greater than one the
*ParFlow* inputs from *pft1* to *pft2* will be repeated the specified *n_cycle* times. (i.e. the total number of time
  steps simulated will be *n_cycle(pft2-pft1)*)
* **ASCII 3-D Particle File Output Frequency (ipwrite):** Controls an ASCII, .3D particle file (not recommended for performance). Refer to the **Outputs** section for the details on this output. If frequency is set to 0 this output will not be written at all, 1 will write the output every timestep, n>1 will write the output every n timesteps (i.e. n=2 writes outputs every other step)
* **VTK Bianary Particle File Output Frequency (ibinpntswrite):** Controls VTK, binary output of particle locations and attributes.  Refer to the **Outputs** section for the details on this output. If frequency is set to 0 this output will not be written at all, 1 will write the output every timestep, n>1 will write the output every n timesteps (i.e. n=2 writes outputs every other step)
* **ASCII ET Output Frequency (etwrite):** Controls ASCII ET output. Refer to the **Outputs** section for the details on this output. If frequency is set to 0 this output will not be written at all, 1 will write the output every timestep, n>1 will write the output every n timesteps (i.e. n=2 writes outputs every other step)
* **VTK Bianary Grid File Output Frequency (icwrite):** Controls VTK, binary grid based outpu, Refer to the **Outputs** section for the details on this output. If frequency is set to 0 this output will not be written at all, 1 will write the output every timestep, n>1 will write the output every n timesteps (i.e. n=2 writes outputs every other step) **##Replace with file names**
* **Velocity Multiplier (V_mult):**  Set to 1.0 for forward particle tracking and -1.0 for backward particle tracking **#Check on backward**
* **CLM Evapotranspiration Flag (clmtrans):** Logical flag (1=True, 0=False) indicating whether evapotranspiration outputs should be read (`pfrunname.out.evapotrans.filenum.pfb`). This is normally generated by *CLM* but can also be generated by *ParFlow* (refer to example **##Fix**).  If this is set to false then no particles are added or removed through recharge and ET. **#CHECK**
* **CLM Variable Read Flag (clmfile):** Logical flag (1=True, 0=False) indicating whether CLM variables should be read to classify recharge events as rain and snow. If false then all recharge events are assumed to be rain. **##Check**
* **Number of Flux Particles (iflux_p_res):** The number of particles to be added per grid cell per ParFLow time step when there is a positive flux into the domain. **#CHECK**
* **Density of Water (denh20):** Density of water [M/L3], used for mass calculations. The units should match the mass and length units set by the *ParFlow* simulation **#CHECK**
* **Molecular Diffusivity(moldiff):** Molecular diffusivity **##Exapand**
* **Numerical Stability Scaler (dtfrac):** Control on the maximum timestep a particle can take to ensure numerical stability (i.e. in one dimension: Particle Timestep= min(pfdt*dtfrac, dtfrac*(dx/Vpx), where Vpx is the velocity in the x direction)  **##Fix**


**Example format for the EcoSLIM input file**
```
SLIM_hillslope   ! SLIM run name, path to ParFlow files follows, then DEM file
"/EcoSLIM/hillslope_clm/hillslope_clm"
"ER_dem.pfb"
20          !nx
5           !ny
5           !nz
20          !particles per cell at start of simulation (-1 = use restart file)
11000000    !np Total
5.0         !dx
0.2         !dy, dz follows
0.1, 0.1, 0.1, 0.1, 0.1
1.0         ! ParFlow DT
1          ! ParFlow t1: ParFlow file number to start from (initial condition is pft1-1)
1752       ! ParFlow t2: ParFlow file number to stop at
1          ! EcoSLIM output start counter 0=pft1
2          ! Time Sequence Repeat [n_cycle*(pft2-pft1)]
0         ! ipwrite frequency, controls an ASCII, .3D particle file not recommended due to poor performance
0         ! ibinpntswrite frequency, controls VTK, binary output of particle locations and attributes
0         !  etwrite frequency, controls ASCII ET output
24        ! icwrite frequency,controls VTK, binary grid based output where particle masses, concentrations
1.0d0       ! velocity multiplier 1.0=forward, -1.0=backward
True         ! CLM Evap Trans Read logical
True           ! CLM Variables Read logical
10          ! number of particles per Evap Trans IC
1000.0      ! density H2O
0.00000414   ! Molecular Diffusivity
0.5d0       ! fraction of Dx/Vx for numerical stability
```

Outputs
-----------------
**Single File Outputs**
1. `runname.out.log`: A log of the settings used for the simulations and any warnings that
occur. **#Check**
2. `runname_particle_restart.bin`: Binary file containing all the particle informaiton for the last time step of the simulation. (116)
3. `runname_exited_particles.bin` (114)
  * Time
  * X
  * Y
  * Z
  * PTIME
  * MASS
  * COMP
  * Exit status
4. `runname_particletrace.3D`: ASCII output writing the particle locaitons at every timestep (X, Y, Z Time) (ipwrite) (Warning- this can be very slow in parallel)  (214)

**Transient File Outputs**
1. `runname_transient_particle.filenum.3D` - ipwrite 1242,
3. `runname_ET_summary.filenum.out.txt` etwrite (1265) - controls ASCII ET output
4. `runname_cgrid_filenum.vtk` icwrite        ! controls VTK, binary grid based output where particle masses, concentrations, ages are mapped to a grid and written every N timesteps
  * 'Concentration'
  * 'Age'
  * 'Mass'
  * 'Comp'
  * Delta'
  * 'ET_Npart'
  * 'ET_Mass'
  * 'ET_Age'
  * 'ET_Comp'787
5. `runname_pnts_filenum.vtk` ibinpnts        !  controls VTK, binary output of particle locations and attributes




Examples and tests contained in the repo
----------------------------------------
The **Examples** folder contains the following test cases. A short description
is provided here. For more details on how to run the examples refer to the
readme files in that directory.
1. **Example name**:Write short description here
2. **Example name**:Write short description here
3. **Example name**:Write short description here
4. **Example name**:Write short description here
