#  This runs the ER shrub PF-CLM hillslope simulation
#  from Maxwell et al EHL 2018
#  the script runs for 5 years with repeated one water-year of
#  hourly forcing to spin up the hillslope
#
# rmaxwell@mines.edu
# May 2018
#
# Import the ParFlow TCL package
#
lappend auto_path $env(PARFLOW_DIR)/bin
package require parflow
namespace import Parflow::*


pfset FileVersion 4

# serial execution
pfset Process.Topology.P 1
pfset Process.Topology.Q 1
pfset Process.Topology.R 1

#---------------------------------------------------------
# Computational Grid
#---------------------------------------------------------
pfset ComputationalGrid.Lower.X           0.0
pfset ComputationalGrid.Lower.Y           0.0
pfset ComputationalGrid.Lower.Z           0.0

pfset ComputationalGrid.NX                20
pfset ComputationalGrid.NY                5
pfset ComputationalGrid.NZ                20

pfset ComputationalGrid.DX	              5.0
pfset ComputationalGrid.DY                0.2
pfset ComputationalGrid.DZ	              0.5

#---------------------------------------------------------
# Domain Geometry
#---------------------------------------------------------
pfset GeomInput.Names                 "domain_input"

#---------------------------------------------------------
# Domain Geometry Input
#---------------------------------------------------------
pfset GeomInput.domain_input.InputType            Box
pfset GeomInput.domain_input.GeomName             domain

#---------------------------------------------------------
# Domain Geometry
#---------------------------------------------------------
pfset Geom.domain.Lower.X                    0.0
pfset Geom.domain.Lower.Y                    0.0
pfset Geom.domain.Lower.Z                    0.0

pfset Geom.domain.Upper.X                  100.0
pfset Geom.domain.Upper.Y                    1.0
pfset Geom.domain.Upper.Z                   10.0

pfset Geom.domain.Patches "x-lower x-upper y-lower y-upper z-lower z-upper"

#--------------------------------------------
# variable dz assignments
#------------------------------------------
pfset Solver.Nonlinear.VariableDz  True
pfset dzScale.GeomNames            domain
pfset dzScale.Type            nzList
pfset dzScale.nzListNumber       20

# 20 layers, starts at 0 for the bottom to  at the top
pfset Cell.0.dzScale.Value 1.0
pfset Cell.1.dzScale.Value 1.0
pfset Cell.2.dzScale.Value 1.0
pfset Cell.3.dzScale.Value 1.0
pfset Cell.4.dzScale.Value 1.0
pfset Cell.5.dzScale.Value 1.0
pfset Cell.6.dzScale.Value 1.0
pfset Cell.7.dzScale.Value 1.0
pfset Cell.8.dzScale.Value 1.0
pfset Cell.9.dzScale.Value 1.0
pfset Cell.10.dzScale.Value 1.0
pfset Cell.11.dzScale.Value 1.0
pfset Cell.12.dzScale.Value 1.0
pfset Cell.13.dzScale.Value 1.0
pfset Cell.14.dzScale.Value 1.0
pfset Cell.15.dzScale.Value 1.0
pfset Cell.16.dzScale.Value 1.0
pfset Cell.17.dzScale.Value 1.0
# 0.5 m * 0.6 = 0.3 m
pfset Cell.18.dzScale.Value .6
# 0.50 m * 0.2 = 0.1m = 10 cm which is default top Noah layer
pfset Cell.19.dzScale.Value 0.2

#-----------------------------------------------------------------------------
# Perm
#-----------------------------------------------------------------------------

pfset Geom.Perm.Names                 "domain"

# Values in m/hour

pfset Geom.domain.Perm.Type            Constant
pfset Geom.domain.Perm.Value           0.05

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
pfset Geom.domain.SpecificStorage.Value 1.0e-5

#-----------------------------------------------------------------------------
# Phases
#-----------------------------------------------------------------------------

pfset Phase.Names "water"

pfset Phase.water.Density.Type	        Constant
pfset Phase.water.Density.Value	        1.0

pfset Phase.water.Viscosity.Type	Constant
pfset Phase.water.Viscosity.Value	1.0

#-----------------------------------------------------------------------------
# Phase sources:
#-----------------------------------------------------------------------------

pfset PhaseSources.water.Type                         Constant
pfset PhaseSources.water.GeomNames                    domain
pfset PhaseSources.water.Geom.domain.Value        0.0


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

# run for 5 years @ 1 hr timesteps, output every timestep
#
pfset TimingInfo.BaseUnit        1.0
pfset TimingInfo.StartCount      0
pfset TimingInfo.StartTime       0.0
pfset TimingInfo.StopTime        43800.0

pfset TimingInfo.DumpInterval    -1
pfset TimeStep.Type              Constant
pfset TimeStep.Value             1.

#-----------------------------------------------------------------------------
# Porosity
#-----------------------------------------------------------------------------

pfset Geom.Porosity.GeomNames          "domain"


pfset Geom.domain.Porosity.Type          Constant
pfset Geom.domain.Porosity.Value         0.2

#-----------------------------------------------------------------------------
# Domain
#-----------------------------------------------------------------------------

pfset Domain.GeomName domain

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------

pfset Phase.RelPerm.Type               VanGenuchten
pfset Phase.RelPerm.GeomNames          "domain"

pfset Geom.domain.RelPerm.Alpha         1.0
pfset Geom.domain.RelPerm.N             2.

#---------------------------------------------------------
# Saturation
#---------------------------------------------------------

pfset Phase.Saturation.Type              VanGenuchten
pfset Phase.Saturation.GeomNames         "domain"

pfset Geom.domain.Saturation.Alpha        1.0
pfset Geom.domain.Saturation.N            2.
pfset Geom.domain.Saturation.SRes         0.2
pfset Geom.domain.Saturation.SSat         1.0


#-----------------------------------------------------------------------------
# Wells
#-----------------------------------------------------------------------------
pfset Wells.Names                           ""

#-----------------------------------------------------------------------------
# Time Cycles
#-----------------------------------------------------------------------------
pfset Cycle.Names "constant"
pfset Cycle.constant.Names           "alltime"
pfset Cycle.constant.alltime.Length  1
pfset Cycle.constant.Repeat         -1


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

pfset Patch.z-upper.BCPressure.Type		      OverlandFlow
pfset Patch.z-upper.BCPressure.Cycle		      "constant"
pfset Patch.z-upper.BCPressure.alltime.Value	      0.0


#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------

pfset TopoSlopesX.Type "Constant"
pfset TopoSlopesX.GeomNames "domain"
pfset TopoSlopesX.Geom.domain.Value 0.1

#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------


pfset TopoSlopesY.Type "Constant"
pfset TopoSlopesY.GeomNames "domain"
pfset TopoSlopesY.Geom.domain.Value 0.00

#---------------------------------------------------------
# Mannings coefficient
#---------------------------------------------------------

pfset Mannings.Type "Constant"
pfset Mannings.GeomNames "domain"
pfset Mannings.Geom.domain.Value 1.e-6

#-----------------------------------------------------------------------------
# Phase sources:
#-----------------------------------------------------------------------------

pfset PhaseSources.Type                         Constant
pfset PhaseSources.GeomNames                    domain
pfset PhaseSources.Geom.domain.Value            0.0

#-----------------------------------------------------------------------------
# Exact solution specification for error calculations
#-----------------------------------------------------------------------------

pfset KnownSolution                                    NoKnownSolution


#-----------------------------------------------------------------------------
# Set solver parameters
#-----------------------------------------------------------------------------

pfset Solver.TerrainFollowingGrid                        True

pfset Solver                                             Richards
pfset Solver.MaxIter                                     2000000

pfset Solver.Nonlinear.MaxIter                           300
pfset Solver.Nonlinear.ResidualTol                       1e-6
pfset Solver.Nonlinear.EtaChoice                         Walker1
pfset Solver.Nonlinear.EtaChoice                         EtaConstant
pfset Solver.Nonlinear.EtaValue                          0.001
pfset Solver.Nonlinear.UseJacobian                       True
pfset Solver.Nonlinear.DerivativeEpsilon                 1e-16
pfset Solver.Nonlinear.StepTol				 1e-20
pfset Solver.Nonlinear.Globalization                     LineSearch
pfset Solver.Linear.KrylovDimension                      20
pfset Solver.Linear.MaxRestart                           2
pfset Solver.Linear.Preconditioner                      PFMG

# no silo output files
pfset Solver.WriteSiloSubsurfData                       False
pfset Solver.WriteSiloPressure                          False
pfset Solver.WriteSiloSaturation                        False
pfset Solver.WriteSiloSlopes                            False
pfset Solver.WriteSiloMask                              False
pfset Solver.WriteSiloEvapTrans                         False
pfset Solver.WriteSiloOverlandSum                       False
pfset Solver.WriteSiloMannings                          False
pfset Solver.WriteSiloSpecificStorage                   False
pfset Solver.WriteSiloCLM                               False
pfset Solver.WriteSiloEvapTrans                         False
pfset Solver.WriteSiloOverlandBCFlux                    False

pfset Solver.PrintOverlandSum                           True
# output velocities in PFB; needed for EcoSLIM
pfset Solver.PrintVelocities                            True
# output EvapTrans matrix in PFB; needed for EcoSLIM input and ET
pfset Solver.PrintEvapTrans                             True
# write all clm output in one file per timestep; EcoSLIM reads this file
# to assign rain / snow
pfset Solver.PrintCLM                                   True
pfset Solver.CLM.SingleFile                             True

# run with CLM
pfset Solver.LSM                                         CLM
pfset Solver.CLM.MetForcing                              1D
pfset Solver.CLM.MetFileName                            ER_forcing.txt
pfset Solver.CLM.MetFilePath                            ../clm_input

# turn off all CLM native binary, log putput
pfset Solver.CLM.Print1dOut                           False
pfset Solver.BinaryOutDir                             False
pfset Solver.WriteCLMBinary                           False
pfset Solver.CLM.WriteLogs                          False

pfset Solver.CLM.CLMDumpInterval                      1

pfset Solver.CLM.EvapBeta                             Linear
pfset Solver.CLM.VegWaterStress                       Saturation
pfset Solver.CLM.ResSat                               0.2
pfset Solver.CLM.WiltingPoint                         0.2
pfset Solver.CLM.FieldCapacity                        1.00
## this key sets the option described in Ferguson, Jefferson, et al ESS 2016
# a setting of 0 (default) will use standard water stress distribution
pfset Solver.CLM.RZWaterStress                           1

# No irrigation
pfset Solver.CLM.IrrigationType                       none


## writing only last daily restarts.  This will be at Midnight GMT and
## starts at timestep 18, then intervals of 24 thereafter
pfset Solver.CLM.WriteLastRST                       True
pfset Solver.CLM.DailyRST                           True

pfset Solver.CLM.RootZoneNZ                         10
pfset Solver.CLM.SoiLayer                           5


#---------------------------------------------------------
# Initial conditions: water pressure
#---------------------------------------------------------

pfset ICPressure.Type                                   HydroStaticPatch
pfset ICPressure.GeomNames                              domain
pfset Geom.domain.ICPressure.RefGeom                    domain
pfset Geom.domain.ICPressure.RefPatch                   z-upper
# hydrostatic IC, WT will be located 2m from ground surface
pfset Geom.domain.ICPressure.Value                      -2.0

#-----------------------------------------------------------------------------
# Run and Unload the ParFlow output files
#-----------------------------------------------------------------------------

file mkdir hillslope_clm_ER_shrub
cd hillslope_clm_ER_shrub

## copy CLM input files
file copy -force ../clm_input/drv_clmin.dat .
file copy -force ../clm_input/drv_vegm.shrub.dat drv_vegm.dat
file copy -force ../clm_input/drv_vegp.dat .

# set name for run; this is prefix for all files output
# and will match EcoSLIM input in slimin.txt
pfrun hillslope_clm_ER_shrub
pfundist hillslope_clm_ER_shrub
cd ..
