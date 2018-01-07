EcoSLIM
=======

**EcoSLIM** private repo for development.  Will be moved to public repo under `/ParFlow` when published

**EcoSLIM** is a Lagrangian, particle-tracking that simulates advective and diffusive movement of water parcels.  This code can be used to simulate age, diagnosing travel times, source water composition and flowpaths.  It integrates seamlessly with **ParFlow-CLM**.

#### Development Team
+ Reed Maxwell <rmaxwell@mines.edu>
+ Mohammad Danesh-Yazdi <danesh@mines.edu>
+ Laura Condon <lecondon@syr.edu>
+ Lindsay Bearup <lbearup@usbr.gov>

To build simply type `make` in the main window, to set number of parallel threads use either
`export OMP_NUM_THREADS=16` for bash or
`setenv OMP_NUM_THREADS 16` for t/c-shell.

To run you will need to have a completed ParFlow simulation and an
EcoSLIM input file that must be named `slimin.txt` and follow the
format described below. Note that this file does not need to be co-located
with the ParFlow simulation.  

To run simply execute `EcoSLIM.exe` from the directory that contains the
`slimin.txt` input file.

`slimin.txt`  Main input file. Includes domain geometry, **ParFlow** timing and input, total number of particles,   initial conditions and information about **CLM**.

### Example format for this file

```
SLIM_hillslope   ! SLIM run name, path to ParFlow files follows
"/EcoSLIM/hillslope_clm/hillslope_clm"
20          !nx
5           !ny
5           !nz
20          !particles per cell at start of simulation
11000000    !np Total
5.0         !dx
0.2         !dy, dz follows
0.1, 0.1, 0.1, 0.1, 0.1
1.0         ! ParFlow DT
1          ! Parflow t1: ParFlow file number to start from (initial condition is pft1-1)
1752       ! Parflow t2: ParFlow file number to stop at
1.0d0       ! velocity multiplier 1.0=forward, -1.0=backward
True        ! CLM Evap Trans
10          ! number of particles per Evap Trans IC
1000.0      ! density H2O
0.00000414   ! Molecular Diffusivity
0.001        ! Fractionation  
0.5d0       ! fraction of Dx/Vx for numerical stability
```
