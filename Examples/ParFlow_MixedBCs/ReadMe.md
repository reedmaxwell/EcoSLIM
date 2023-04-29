ParFlow MixedBC Example
-------------------------
This is a 2D hillslope example with fixed head boundaries on the sides and constant flux across the top of the domain. It demonstrates use of the ecoslim *velfile* option, but note that this example must be run with *ParFlow* >= 3.10.0

First the *ParFlow* domain is run to steady state then a transient simulation is completed which is used by *EcoSLIM*.

Because the same forcing is used for the transient simulation as the steady state simulation in this example, the same result could have been achieved using only the steady state *ParFlow* outputs as input to *EcoSLIM* and increasing *n_cycle* in `slimin.txt` (i.e. repeating the steady state outputs for *EcoSLIM* as many times as desired). Here we complete a transient simulation to demonstrate the `slimin.txt` setup for transient *ParFlow* simulations. The example is currently set to run over the entire timeseries of *ParFlow* outputs but the users can experiment with changing *pft1*, *pft2* and *n_cycle* in `slimin.txt`.

## Steps to run the example
1. Copy this directory from the GitHub repository to where you want to work

2. Run the spinup tcl script: `tclsh mixedBCs_hillslope2D_spinup.tcl`
  * This will create the directory `mixedBCs_hillslope2D_spinup` that should have 20 timesteps of output.

3. Run the transient tcl script: `tclsh mixedBCs_hillslope2D_transient.tcl`
  * This will create the directory 'mixedBCs_hillslope2D_transient'
  * It will copy the last time step from the spinup run to this directory and use it as the initial condition for the transient run
  * When the run finishes it will create the 'EcoSLIM_mixedBCs_hillslope2D' folder and copy the EcoSLIM input into it

4. Run EcoSLIM:
  * Change directory to where the *EcoSLIM* input file sits `cd EcoSLIM_mixedBCs_hillslope2D`
  * Run the `EcoSLIM.exe` executable from this directory `.\path_to_EcoSLIM\EcoSLIM.exe`
  * This should give you 1000 timesteps of VTK outputs as well as the gridded
 text outputs.
