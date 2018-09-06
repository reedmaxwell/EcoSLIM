This directory contain hillslope simulation examples that demonstrate EcoSLIM

** Paper Cases : these are the cases from [Maxwell et al Ecohydrology](https://doi.org/10.1002/eco.2042)

 These cases run in two parts, first the ParFlow tcl script, then the EcoSLIM run.  They correspond to the *LW* and *ER* forcing and the *Shrub* and *Trees* cases.  The ParFlow simulations run for five years to spin up, the EcoSLIM cases for 20 years, forced by only the last year of PF-CLM output.

 first run a ParFlow simulation (e.g.):
```
tclsh hillslope_clm.ER_shrub.tcl
```
then cd into the corresponding directory (e.g.):
```
cd SLIM_ER_hillslope_20y_shrub
```
set threads (in this case 4):
```
export OMP_NUM_THREADS=4
```

run EcoSLIM
```
../../../../EcoSLIM.exe
```
