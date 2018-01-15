
lappend   auto_path $env(PARFLOW_DIR)/bin
package   require parflow
namespace import Parflow::*

set       vx         [pfload -pfb ./spinup0.00125/spinup0.00125.out.velx.00020.pfb]
#pfsetgrid {209 268 50} {534825.340 4584570.246 0.0} {812.904 812.904 1.0} $DZICP

set       vy         [pfload -pfb ./spinup0.00125/spinup0.00125.out.vely.00020.pfb]
#pfsetgrid {209 268 50} {534825.340 4584570.246 0.0} {812.904 812.904 1.0} $DZICP

set       vz         [pfload -pfb ./spinup0.00125/spinup0.00125.out.velz.00020.pfb]
#pfsetgrid {209 268 50} {534825.340 4584570.246 0.0} {812.904 812.904 1.0} $DZICP


# Silo...
pfsave $vx -silo spinup0.00125.out.velx.00020.silo
pfsave $vy -silo spinup0.00125.out.vely.00020.silo
pfsave $vz -silo spinup0.00125.out.velz.00020.silo
