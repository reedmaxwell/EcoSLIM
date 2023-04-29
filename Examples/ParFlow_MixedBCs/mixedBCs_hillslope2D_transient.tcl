# Mixed boundary condition test case for EcoSLIM
# Demonstrates incorporation of ParFlow DirEquilRefPatch and 
# FluxConst boundary conditions across a hillslope

set tcl_precision 17

# Get the run name from the file name
set file_name [info script]
set end_index [expr [string length $file_name] - 5]
set run_name [string range $file_name 0 $end_index]
file mkdir $run_name
cd $run_name

# Get ic_file
set end_index [expr [string length $run_name] - 11]
set run_name_prefix [string range $run_name 0 $end_index]
set ic_file [format %s_spinup.out.press.00020.pfb $run_name_prefix] 

#
# Import the ParFlow TCL package
#
lappend auto_path $env(PARFLOW_DIR)/bin
package require parflow
namespace import Parflow::*

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
pfset ComputationalGrid.NY                1
pfset ComputationalGrid.NZ                40

pfset ComputationalGrid.DX	             1.0
pfset ComputationalGrid.DY               1.0
pfset ComputationalGrid.DZ	             .1

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
pfset Geom.domain.Upper.Y                        1.0
pfset Geom.domain.Upper.Z                        4.0
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
pfset TimingInfo.StopTime        1000
pfset TimingInfo.DumpInterval    -1

pfset TimeStep.Type              Constant
pfset TimeStep.Value    1.0

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

pfset Patch.x-lower.BCPressure.Type		      DirEquilRefPatch
pfset Patch.x-lower.BCPressure.Cycle		      "constant"
pfset Patch.x-lower.BCPressure.RefGeom	      domain
pfset Patch.x-lower.BCPressure.RefPatch	      z-lower
pfset Patch.x-lower.BCPressure.alltime.Value	      1.5

pfset Patch.y-lower.BCPressure.Type		      FluxConst
pfset Patch.y-lower.BCPressure.Cycle		      "constant"
pfset Patch.y-lower.BCPressure.alltime.Value	      0.0

pfset Patch.z-lower.BCPressure.Type		      FluxConst
pfset Patch.z-lower.BCPressure.Cycle		      "constant"
pfset Patch.z-lower.BCPressure.alltime.Value	      0.0

pfset Patch.x-upper.BCPressure.Type		      DirEquilRefPatch
pfset Patch.x-upper.BCPressure.Cycle		      "constant"
pfset Patch.x-upper.BCPressure.RefGeom	      domain
pfset Patch.x-upper.BCPressure.RefPatch	      z-lower
pfset Patch.x-upper.BCPressure.alltime.Value	      3.0

pfset Patch.y-upper.BCPressure.Type		      FluxConst
pfset Patch.y-upper.BCPressure.Cycle		      "constant"
pfset Patch.y-upper.BCPressure.alltime.Value	      0.0

pfset Patch.z-upper.BCPressure.Type		      FluxConst
pfset Patch.z-upper.BCPressure.Cycle		      "constant"
pfset Patch.z-upper.BCPressure.alltime.Value	      -0.0005

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
pfset Solver.MaxIter                                     25000
pfset OverlandFlowDiffusive  0


pfset Solver.Nonlinear.MaxIter                           20
pfset Solver.Nonlinear.ResidualTol                       1e-7
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
pfset Solver.Linear.Preconditioner.PCMatrixType     FullJacobian
pfset Solver.PrintSubsurf				False
pfset  Solver.Drop                                      1E-20
pfset Solver.AbsTol                                     1E-7

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
set spinup_dir [format %s_spinup $run_name_prefix]
file copy -force ../$spinup_dir/$ic_file "./ICPress.Spinup.pfb"

# set water table to be at the bottom of the domain, the top layer is initially dry
pfset ICPressure.Type                                   PFBFile
pfset ICPressure.GeomNames                              domain
pfset Geom.domain.ICPressure.FileName                   ICPress.Spinup.pfb


#-----------------------------------------------------------------------------
# Run and Unload the ParFlow output files
#-----------------------------------------------------------------------------


pfdist ICPress.Spinup.pfb
pfrun $run_name
pfundist $run_name
pfundist ICPress.Spinup.pfb

cd ..

file mkdir EcoSLIM_$run_name_prefix
cd EcoSLIM_$run_name_prefix
file copy -force ../mixedBCs_slimin.txt slimin.txt

