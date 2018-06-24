
set tcl_precision 17

set runname flux_hillslope
#
# Import the ParFlow TCL package
#
lappend auto_path $env(PARFLOW_DIR)/bin 
package require parflow
namespace import Parflow::*



set fileId [open dem.sa w]
puts $fileId "20 5 1"
for { set jj 0 } { $jj < 5 } { incr jj } {
for { set ii 0 } { $ii < 20 } { incr ii } {
puts $fileId [expr ((($ii+1)*5.0-2.5)*0.1)] 
}
}


close $fileId

set       dem         [pfload -sa dem.sa]
pfsetgrid {20 5 1} {0.0 0.0 0.0} {5.0 0.2 0.05} $dem
pfsave $dem -pfb hillslope.dem.pfb 
