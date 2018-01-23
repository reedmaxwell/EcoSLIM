#!/bin/bash

#Run ParFlow for rage of start dates
#for fluxrate in 0.00125 0.006 0.0125 0.0250 0.0375 0.050
#do
#        echo "$d"
#        #tclsh steadyflux_hillslope2D_spinup.tcl $fluxrate
#        mkdir "EcoSLIM$fluxrate"
#        cp
#done

# se
name=Mytest
Rscript WriteSlimin.R $name
