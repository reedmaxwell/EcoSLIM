# overland flux test case for EcoSLIM
# provides some rainfall and then ET
# over a hillslope

set tcl_precision 17

set runname flux_hillslope
#
# Import the ParFlow TCL package
#
lappend auto_path $env(PARFLOW_DIR)/bin 
package require parflow
namespace import Parflow::*


file mkdir flux_hillslope
cd flux_hillslope

pfset FileVersion 4

pfset Process.Topology.P        1
pfset Process.Topology.Q        1
pfset Process.Topology.R        1

#---------------------------------------------------------
# Computational Grid
#---------------------------------------------------------
pfset ComputationalGrid.Lower.X           0.0
pfset ComputationalGrid.Lower.Y           0.0
pfset ComputationalGrid.Lower.Z           0.0

pfset ComputationalGrid.NX                30
pfset ComputationalGrid.NY                3
pfset ComputationalGrid.NZ                30

pfset ComputationalGrid.DX	         1.0
pfset ComputationalGrid.DY               1.0
pfset ComputationalGrid.DZ	            .05

#---------------------------------------------------------
# The Names of the GeomInputs
#---------------------------------------------------------
pfset GeomInput.Names                 "domaininput"

pfset GeomInput.domaininput.GeomName  domain
pfset GeomInput.domaininput.InputType  Box 

#---------------------------------------------------------
# Domain Geometry 
#---------------------------------------------------------
pfset Geom.domain.Lower.X                        0.0
pfset Geom.domain.Lower.Y                        0.0
pfset Geom.domain.Lower.Z                        0.0
 
pfset Geom.domain.Upper.X                        30.0
pfset Geom.domain.Upper.Y                        3.0
pfset Geom.domain.Upper.Z                          1.5
pfset Geom.domain.Patches             "x-lower x-upper y-lower y-upper z-lower z-upper"


#-----------------------------------------------------------------------------
# Perm
#-----------------------------------------------------------------------------

pfset Geom.Perm.Names                 "domain"

# Values in m/hour#

pfset Geom.domain.Perm.Type            Constant
pfset Geom.domain.Perm.Value           0.1


pfset Perm.TensorType               TensorByGeom

pfset Geom.Perm.TensorByGeom.Names  "domain"

pfset Geom.domain.Perm.TensorValX  1.0d0
pfset Geom.domain.Perm.TensorValY  1.0d0
pfset Geom.domain.Perm.TensorValZ  1.0d0

#-----------------------------------------------------------------------------
# Specific Storage
#-----------------------------------------------------------------------------

pfset SpecificStorage.Type            Constant
pfset SpecificStorage.GeomNames       "domain"
pfset Geom.domain.SpecificStorage.Value 1.0e-4

#-----------------------------------------------------------------------------
# Phases
#-----------------------------------------------------------------------------

pfset Phase.Names "water"

pfset Phase.water.Density.Type	        Constant
pfset Phase.water.Density.Value	        1.0

pfset Phase.water.Viscosity.Type	Constant
pfset Phase.water.Viscosity.Value	1.0

#-----------------------------------------------------------------------------
# Contaminants
#-----------------------------------------------------------------------------

pfset Contaminants.Names			""

#-----------------------------------------------------------------------------
# Retardation
#-----------------------------------------------------------------------------

pfset Geom.Retardation.GeomNames           ""

#-----------------------------------------------------------------------------
# Gravity
#-----------------------------------------------------------------------------

pfset Gravity				1.0

#-----------------------------------------------------------------------------
# Setup timing info
#-----------------------------------------------------------------------------

# 
pfset TimingInfo.BaseUnit        0.1
pfset TimingInfo.StartCount      0
pfset TimingInfo.StartTime       0.0
pfset TimingInfo.StopTime        2.5
pfset TimingInfo.DumpInterval    -1
pfset TimeStep.Type              Constant
pfset TimeStep.Value             0.1
 
#-----------------------------------------------------------------------------
# Porosity
#-----------------------------------------------------------------------------

pfset Geom.Porosity.GeomNames          "domain"

pfset Geom.domain.Porosity.Type          Constant
pfset Geom.domain.Porosity.Value         0.25


#-----------------------------------------------------------------------------
# Domain
#-----------------------------------------------------------------------------

pfset Domain.GeomName domain

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------

pfset Phase.RelPerm.Type               VanGenuchten
pfset Phase.RelPerm.GeomNames          "domain"

pfset Geom.domain.RelPerm.Alpha         6.0
pfset Geom.domain.RelPerm.N             2. 

#---------------------------------------------------------
# Saturation
#---------------------------------------------------------

pfset Phase.Saturation.Type              VanGenuchten
pfset Phase.Saturation.GeomNames         "domain"

pfset Geom.domain.Saturation.Alpha        6.0
pfset Geom.domain.Saturation.N            2.
pfset Geom.domain.Saturation.SRes         0.02
pfset Geom.domain.Saturation.SSat         1.0



#-----------------------------------------------------------------------------
# Wells
#-----------------------------------------------------------------------------
pfset Wells.Names                           ""

#-----------------------------------------------------------------------------
# Time Cycles
#-----------------------------------------------------------------------------
pfset Cycle.Names "constant"
pfset Cycle.constant.Names              "alltime"
pfset Cycle.constant.alltime.Length      1
pfset Cycle.constant.Repeat             -1

 
#-----------------------------------------------------------------------------
# Boundary Conditions: Pressure
#-----------------------------------------------------------------------------
pfset BCPressure.PatchNames                   [pfget Geom.domain.Patches]

pfset Patch.x-lower.BCPressure.Type		      FluxConst
pfset Patch.x-lower.BCPressure.Cycle		      "constant"
pfset Patch.x-lower.BCPressure.alltime.Value	      0.0

pfset Patch.y-lower.BCPressure.Type		      FluxConst
pfset Patch.y-lower.BCPressure.Cycle		      "constant"
pfset Patch.y-lower.BCPressure.alltime.Value	      0.0

pfset Patch.z-lower.BCPressure.Type		      FluxConst
pfset Patch.z-lower.BCPressure.Cycle		      "constant"
pfset Patch.z-lower.BCPressure.alltime.Value	      0.0

pfset Patch.x-upper.BCPressure.Type		      FluxConst
pfset Patch.x-upper.BCPressure.Cycle		      "constant"
pfset Patch.x-upper.BCPressure.alltime.Value	      0.0

pfset Patch.y-upper.BCPressure.Type		      FluxConst
pfset Patch.y-upper.BCPressure.Cycle		      "constant"
pfset Patch.y-upper.BCPressure.alltime.Value	      0.0

## overland flow boundary condition with very heavy rainfall then slight ET
pfset Patch.z-upper.BCPressure.Type		      OverlandFlow
pfset Patch.z-upper.BCPressure.Cycle		      "constant"
pfset Patch.z-upper.BCPressure.alltime.Value	      0.0

pfset Solver.EvapTransFileTransient True
pfset Solver.EvapTrans.FileName Forcing 

set fileId [open flux.rain.sa w]
puts $fileId "30 3 30"
for { set kk 0 } { $kk < 30 } { incr kk } {
for { set jj 0 } { $jj < 3 } { incr jj } {
for { set ii 0 } { $ii < 30 } { incr ii } {

	if {$kk == 29} {
		# convert from 0.05 m / h of rainfall to flux 1/t units; note negative positive flip
		# and divide by dz
		puts $fileId "0.0" 
		#puts $fileId  [expr (0.05/0.05)] 
	} else {
		puts $fileId "0.0"  }
}
}
}

close $fileId

set fileId [open flux.rec.sa w]
puts $fileId "30 3 30"
for { set kk 0 } { $kk < 30 } { incr kk } {
for { set jj 0 } { $jj < 3 } { incr jj } {
for { set ii 0 } { $ii < 30 } { incr ii } {

	# top 5 cells have ET like a simple root zone
		if {$kk == 26} {
		# convert from 0.001 m / h of et to flux 1/t units; note negative positive flip
		# and divide by dz
		puts $fileId [expr (-0.005/0.05)] 
	} else {
		puts $fileId "0.0"  }

	}
}
}

close $fileId

set       flux         [pfload -sa flux.rain.sa]
pfsetgrid {30 3 30} {0.0 0.0 0.0} {10.0 10.0 0.05} $flux
pfsave $flux -pfb Forcing.00000.pfb 
pfsave $flux -silo Forcing.00000.silo 
pfsave $flux -pfb Forcing.00001.pfb 
pfsave $flux -pfb Forcing.00002.pfb 
pfsave $flux -pfb Forcing.00003.pfb 
pfsave $flux -pfb Forcing.00004.pfb
pfsave $flux -pfb Forcing.00005.pfb 
pfsave $flux -pfb Forcing.00006.pfb 
pfsave $flux -pfb Forcing.00007.pfb 
pfsave $flux -pfb Forcing.00008.pfb 
pfsave $flux -pfb Forcing.00009.pfb 
pfsave $flux -pfb Forcing.00010.pfb 
set       flux         [pfload -sa flux.rec.sa]
pfsetgrid {30 3 30} {0.0 0.0 0.0} {10.0 10.0 0.05} $flux
pfsave $flux -pfb Forcing.00011.pfb 
pfsave $flux -pfb Forcing.00012.pfb 
pfsave $flux -pfb Forcing.00013.pfb 
pfsave $flux -pfb Forcing.00014.pfb 
pfsave $flux -pfb Forcing.00015.pfb 
pfsave $flux -pfb Forcing.00016.pfb 
pfsave $flux -pfb Forcing.00017.pfb 
pfsave $flux -pfb Forcing.00018.pfb 
pfsave $flux -pfb Forcing.00019.pfb 
pfsave $flux -pfb Forcing.00020.pfb 
pfsave $flux -pfb Forcing.00021.pfb 
pfsave $flux -pfb Forcing.00022.pfb 
pfsave $flux -pfb Forcing.00023.pfb 
pfsave $flux -pfb Forcing.00024.pfb 
pfsave $flux -pfb Forcing.00025.pfb 

pfdist  Forcing.00000.pfb 
pfdist  Forcing.00001.pfb 
pfdist  Forcing.00002.pfb 
pfdist  Forcing.00003.pfb 
pfdist  Forcing.00004.pfb 
pfdist  Forcing.00005.pfb 
pfdist  Forcing.00006.pfb 
pfdist  Forcing.00007.pfb 
pfdist  Forcing.00008.pfb 
pfdist  Forcing.00009.pfb 
pfdist  Forcing.00010.pfb 
pfdist  Forcing.00011.pfb 
pfdist  Forcing.00012.pfb 
pfdist  Forcing.00013.pfb 
pfdist  Forcing.00014.pfb 
pfdist  Forcing.00015.pfb 
pfdist  Forcing.00016.pfb 
pfdist  Forcing.00017.pfb 
pfdist  Forcing.00018.pfb 
pfdist  Forcing.00019.pfb 
pfdist  Forcing.00020.pfb 
pfdist  Forcing.00021.pfb 
pfdist  Forcing.00022.pfb 
pfdist  Forcing.00023.pfb 
pfdist  Forcing.00024.pfb 
pfdist  Forcing.00025.pfb 



#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------

pfset TopoSlopesX.Type "Constant"
pfset TopoSlopesX.GeomNames "domain"
pfset TopoSlopesX.Geom.domain.Value 0.00

#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------


pfset TopoSlopesY.Type "Constant"
pfset TopoSlopesY.GeomNames "domain"
pfset TopoSlopesY.Geom.domain.Value 0.0

#---------------------------------------------------------
# Mannings coefficient 
#---------------------------------------------------------

pfset Mannings.Type "Constant"
pfset Mannings.GeomNames "domain"
pfset Mannings.Geom.domain.Value 1.e-6

#-----------------------------------------------------------------------------
# Phase sources:
#-----------------------------------------------------------------------------

pfset PhaseSources.water.Type                         Constant
pfset PhaseSources.water.GeomNames                    domain
pfset PhaseSources.water.Geom.domain.Value        0.0

#-----------------------------------------------------------------------------
# Exact solution specification for error calculations
#-----------------------------------------------------------------------------

pfset KnownSolution                                    NoKnownSolution


#-----------------------------------------------------------------------------
# Set solver parameters
#-----------------------------------------------------------------------------

pfset Solver                                             Richards
pfset Solver.MaxIter                                     2500
pfset OverlandFlowDiffusive  0


pfset Solver.Nonlinear.MaxIter                           20
pfset Solver.Nonlinear.ResidualTol                       1e-9
pfset Solver.Nonlinear.EtaChoice                         EtaConstant
pfset Solver.Nonlinear.EtaValue                          0.01
pfset Solver.Nonlinear.UseJacobian                       False
pfset Solver.Nonlinear.UseJacobian                      True 
pfset Solver.Nonlinear.DerivativeEpsilon                 1e-8
pfset Solver.Nonlinear.StepTol				 1e-20
pfset Solver.Nonlinear.Globalization                     LineSearch
pfset Solver.Linear.KrylovDimension                      20
pfset Solver.Linear.MaxRestart                           2

pfset Solver.Linear.Preconditioner                       PFMG
#pfset Solver.Linear.Preconditioner                       MGSemi 
pfset Solver.Linear.Preconditioner.PCMatrixType     FullJacobian
pfset Solver.PrintSubsurf				False
pfset  Solver.Drop                                      1E-20
pfset Solver.AbsTol                                     1E-9
 
pfset Solver.WriteSiloSubsurfData True
pfset Solver.WriteSiloPressure True
pfset Solver.WriteSiloSaturation True
pfset Solver.WriteSiloEvapTrans                         True


######
## Make sure we write PFB output for EcoSLIM
#
pfset Solver.PrintVelocities   			    True
pfset Solver.PrintEvapTrans                         True



#---------------------------------------------------------
# Initial conditions: water pressure
#---------------------------------------------------------

# set water table to be at the bottom of the domain, the top layer is initially dry
pfset ICPressure.Type                                   HydroStaticPatch
pfset ICPressure.GeomNames                              domain
pfset Geom.domain.ICPressure.Value                      -1.0

pfset Geom.domain.ICPressure.RefGeom                    domain
pfset Geom.domain.ICPressure.RefPatch                   z-upper

#-----------------------------------------------------------------------------
# Run and Unload the ParFlow output files
#-----------------------------------------------------------------------------


pfrun $runname
pfundist $runname

cd ..


