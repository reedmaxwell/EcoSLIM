ParFlow BoolSourceFile Example
-------------------------
This is a 2D hillslope example with fixed head boundaries on the sides and constant flux across the top of the domain. It demonstrates use of the ecoslim *velfile* and *Boolfile* options, but note that itmust be run with *ParFlow* >= 3.10.0

First the *ParFlow* domain is run to steady state then a transient simulation is completed which is used by *EcoSLIM*. 

Within *EcoSLIM*, initial particle locations are limited to a portion of the domain using boolean particle source files. Initially, particles are located in a rectangle with 15<=X<=26 m and -3.1<=Z-1.0 m. For all subsequent time steps, particles are only added to the upper half of the domain (X >= 10 m).

## Steps to run the example
1. Copy this directory from the GitHub repository to where you want to work

2. Run the spinup tcl script: `tclsh BoolSourceFile_hillslope2D_spinup.tcl`
  * This will create the directory `BoolSourceFile_hillslope2D_spinup` that should will have 20 timesteps of output.

3. Run the transient tcl script: `tclsh BoolSourceFile_hillslope2D_transient.tcl`
  * This will create the directory `BoolSourceFile_hillslope2D_transient`
  * It will copy the last time step from the spinup run to this directory and use it as the initial condition for the transient run
  * When the run finishes it will create the `EcoSLIM_BoolSourceFile_hillslope2D` folder and copy the EcoSLIM input into it
  * It will also create boolean particle source files for each time step (`PSourceBool.#####.pfb`)

4. Run EcoSLIM:
  * Change directory to where the *EcoSLIM* input file sits `cd EcoSLIM_BoolSourceFile_hillslope2D`
  * Run the `EcoSLIM.exe` executable from this directory `.\path_to_EcoSLIM\EcoSLIM.exe`
  * This should give you 1000 timesteps of VTK outputs as well as the gridded
 text outputs.
