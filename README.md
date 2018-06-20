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

*Maxwell, R.M., L.E. Condon, M. Danesh-Yazdi and L.A. Bearup, *Exploring source water mixing and transient residence time distributions of outflow and  evapotranspiration with an integrated hydrologic model and Lagrangian particle tracking approach*, in review, 2018.*

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

Required Inputs
--------------------
**Mandatory Inputs**
1. `slimin.txt`: EcoSLIM input file (See following section for description)
2.

**Optional Inputs**
1. DEM for ParFlow simulation: If a DEM is provided... **##FILL in DEM Details**
2.

EcoSLIM Parameters
--------------------
The `slimin.txt` file contains all of the settings for an EcoSLIM Simulation.
An example of this file is provided below and other examples are also provided
in the **Examples** folder. Inputs *must* be provided in the oder they appear in the
template shown here.

Here we describe each parameter and what they do:
1. EcoSLIM Runname:
⋅⋅⋅The runname for all of the file outputs
* ParFlow Runname: The runname used for the ParFlow simulations. If the **ParFlow**
outputs are not located in the same directory as the **EcoSLIM** run this should also
include the directory path to the **ParFlow** simulation as shown in the example below.
* DEM File Name: The name of the DEM file for the **ParFlow** simulations. If you want
to provide a dem you can leave this line blank **##CLARIFY THIS**
* **ParFlow** nx: Number of grid cells in the x-direction for ParFlow domain
* **ParFlow** ny:  Number of grid cells in the y-direction for ParFlow domain
* **ParFlow** nz:  Number of grid cells in the z-direction for ParFlow domain
* Initial Particles: If a positive integer is provided this will be the number of
particles placed in every grid cell at the start of the simulation. To start from a previous
**EcoSLIM** output set this value to -1. In this case particles will be initialized from the
`runname_particle_restart.bin` file (refer to the outputs section for details on this file).  **##CLARIFY how this gets overwritten**
* Total Particles: The total number of particles allowed in the simulation. If the particle count
exceeds this at any point (i.e. through particle addition with initial conditions or rainfall events)
the simulation will exit.
* **ParFlow** dx: **ParFlow** grid cell size in the x-direction
* **ParFlow** dy: **ParFlow** grid cell size in the y-direction
* **ParFlow** dz: **ParFlow** grid cell size in the x-direction. This should be a list separated
by comas that is nz long (refer to example below).
* **ParFlow** dt: The time step used for the ParFlow simulation. Note that currently EcoSLIM assumes
a constant time step so this is not compatible with growth times steps in **ParFlow**
* Initial ParFlow Timestep:



including/defining:
* Number of particles per cell at start of simulation
* Total number of particles to be tracked
* ParFlow timestep
* ParFlow file number to start from
* ParFlow file number to stop at
* Velocity multiplier: 1.0=forward tracking, -1.0=backward tracking
* A logical parameter that turns CLM Evap Trans reading on / off
* Number of particles entering domain via Evap Trans
* Density of water  
* A logical parameter that turns CLM output reading on / off
* Molecular Diffusivity
* Fraction of Dx/Vx (also Dy/Vy and Dz/Vz) for numerical stability
* Number of concentration constituents


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

Model Outputs
-----------------
1. `runname.out.log`: A log of the settings used for the simulations and any warnings that
occured **CONFIRM**
2. `runname_particle_restart.bin`



Examples and tests contained in the repo
----------------------------------------
The **Examples** folder contains the following test cases. A short description
is provided here. For more details on how to run the examples refer to the
readme files in that directory.
1. **Example name**:Write short description here
2. **Example name**:Write short description here
3. **Example name**:Write short description here
4. **Example name**:Write short description here
