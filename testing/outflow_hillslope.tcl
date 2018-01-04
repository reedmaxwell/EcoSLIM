#
# Import the ParFlow TCL package
#
lappend auto_path $env(PARFLOW_DIR)/bin
package require parflow
namespace import Parflow::*

pfset FileVersion 4
set tclprecision 12

## open a file
set filenum [open hillslope_clm_out.txt w]

cd hillslope_clm
set runname hillslope_clm

puts $filenum "Time(hr) Runoff(m^3/hr) P-ET(m^3/hr) Well1(m) Well2(m) Well3(m)"
for {set i 1} {$i <= 1752} {incr i} {

## first get total outflow
	set filename [format "%s.out.overlandsum.%05d.silo" $runname $i]
	set surface_runoff2 [pfload $filename]
	set total_surface_runoff2 [pfsum $surface_runoff2]
## get P-ET
	set filename [format "%s.out.evaptrans.%05d.silo" $runname $i]
	set PmET [pfload $filename]
	set total_PmET [pfsum $PmET]

## then get well OBS
	set  filename [format "%s.out.press.%05d.silo" $runname $i]
	set pressure [pfload $filename]
        set well1 [pfgetelt $pressure 1 0 10]
        set well2 [pfgetelt $pressure 5 0 10]
        set well3 [pfgetelt $pressure 10 0 10]
## write to a file, make it look neater in columns
	set outformat "%8i %5e %- 5e %5e %5e %5e"
	puts $filenum [format $outformat $i $total_surface_runoff2 $total_PmET $well1 $well2 $well3]
}

close $filenum

cd ..
