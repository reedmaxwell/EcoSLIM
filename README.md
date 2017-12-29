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

`slimin.txt`  Main input file. Includes domain geometry, **ParFlow** timing and input, total number of particles,   initial conditions and information about **CLM**.
