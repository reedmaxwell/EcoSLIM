#
# executable & I/O
#

set slim_executable $env(SLIM_DIR)/bin/SLIM.exe

#
#
# Write to the Parameter file:
#
set fileId [open $parameter_file w 0600]

puts $fileId "\"$run_name\"						! Title of Run"
puts $fileId "\"$log_file\"						! Log File Name"

if {$saturated == "yes"} {
puts $fileId "1"
} else {
puts $fileId "0"}


if {$v_type == "vbin"} {
puts $fileId "1							! reading V's from vbin type file"
puts $fileId "\"$vel_file\"					! Vbin File Name" }


if {$v_type == "calc"} {
puts $fileId "2							! calcing V's internally, steady state"
puts $fileId "\"$kx_file\"								! Kx File Name" 
puts $fileId "\"$ky_file\"								! Ky File Name" 
puts $fileId "\"$kz_file\"								! Kz File Name" 
## @RMM 8-6-08 added to make consistent w/ varsat
if {$saturated != "yes"} {
puts $fileId "\"$vga_file\"								! VG alpha File Name" 
puts $fileId "\"$vgn_file\"								! VG n File Name" 
puts $fileId "\"$sres_file\"								! VG Sres File Name"
} 
puts $fileId "$press								! head (0) or pressure (1) flag" 
puts $fileId "\"$head_file\"								! Head File Name" 
if {$phi_type == "constant"} {
puts $fileId "1										! const phi"
puts $fileId "$phi									! phi value \[-\]"
}
if {$phi_type == "PFBFile"} {
puts $fileId "2										! const phi"
puts $fileId "\"$phi_file\"									! phi value \[-\]"
}
}
if {$v_type == "calc_trans"} {
puts $fileId "3							! calcing v's internally, transient"
puts $fileId "\"$kx_file\"								! Kx File Name" 
puts $fileId "\"$ky_file\"								! Ky File Name" 
puts $fileId "\"$kz_file\"								! Kz File Name" 
if {$saturated != "yes"} {
puts $fileId "\"$vga_file\"								! VG alpha File Name" 
puts $fileId "\"$vgn_file\"								! VG n File Name" 
puts $fileId "\"$sres_file\"								! VG Sres File Name" }
puts $fileId "$press								! head (0) or pressure (1) flag" 
puts $fileId "\"$time_file\"								! File w/ timesteps " 
puts $fileId "\"$head_list_file\"								! File w/  pressure or head" 
if {$phi_type == "constant"} {
puts $fileId "1										! const phi"
puts $fileId "$phi									! phi value \[-\]"
}
if {$phi_type == "PFBFile"} {
puts $fileId "2										! const phi"
puts $fileId "\"$phi_file\"									! phi value \[-\]"
}
}
puts $fileId "$nx 						! number of x-nodes in the domain \[-\]"
puts $fileId "$ny 						! number of y-nodes in the domain \[-\]"
puts $fileId "$nz 						! number of z-nodes in the domain \[-\]"
puts $fileId "$dx 						! delta-x \[m\]"
puts $fileId "$dy 						! delta-y \[m\]"
puts $fileId "$dz 						! delta-z \[m\]"
puts $fileId "$num_constituents 				! number of constituents"
for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {
puts $fileId "$half_life($jj)				        ! radioactive half life,  \[d\]" 
}
puts $fileId "$alpha_l						! longitudinal disersivity, alpha_l, \[m\]"
puts $fileId "$alpha_t						! transverse disersivity, alpha_t, \[m\]"

for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {
if {$sorption_type($jj) == "constant"} {
puts $fileId "$sorption_value($jj)					! Linear Retardation Coeff, R, \[-\]" 
}
if {$sorption_type($jj) == "min_file"} {
puts $fileId "-9999					! min_file" 
puts $fileId "\"$mineralization_file($jj)\""
}
if {$sorption_type($jj) == "kd_file"} {
puts $fileId "-1					! kd_file" 
puts $fileId "\"$kd_file($jj)\""
puts $fileId "$bulk_density($jj) 				! bulk density"

}
if {$sorption_type($jj) == "R_file"} {
puts $fileId "-2					! R_file" 
puts $fileId "\"$R_file($jj)\""

}

}

if {$fast_kin == "yes"} {
puts $fileId "1"
} else {
puts $fileId "0"}


for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {
if {$attach_type($jj) == "constant"} {
puts $fileId "$attachment_value($jj)					! Attachment Coeff, Katt, \[1/d\]" 
puts $fileId "$detachment_value($jj)					! Detachment Coeff, Kdet, \[1/d\]" 
}
if {$attach_type($jj) == "file"} {
puts $fileId "-9999					! file" 
puts $fileId "-9999					! file" 
puts $fileId "\"$att_det_file($jj)\""
puts $fileId "\"$perm_catagory_file($jj)\""
}

#if {$detach_type($jj) == "constant"} {
#puts $fileId "$detachment_value($jj)					! Detachment Coeff, Kdet, \[1/d\]" 
#}
#if {$detach_type($jj) == "file"} {
#puts $fileId "-9999					! file" 
#}

}


puts $fileId "$time_increment					! time increment until reporting particle vales, Tnext \[d\]"
puts $fileId "$num_increments					! number of time increments, nt, (total runtime=nt*tnext in days)"

if {$write_well_breakthrough == "yes"} {if {$well_overwrite == "no"} {puts $fileId "1				! appending well BTC to file:"
} else { puts $fileId "2				! writing out well BTC to file:"}
} else { puts $fileId "0				! not writing out well BTC"}

for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {
puts $fileId "\"$well_out_file($jj)\""
}

puts $fileId "$well_breakthrough_dt				! time discretization for wells \[days\]"
puts $fileId "$number_well_steps				! number of time discretizations for wells"
if {$write_concentrations == "yes"} {puts $fileId "1				! writing out concentrations using header:"
for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {
puts $fileId "\"$concentration_header($jj)\""
}
} elseif {$write_concentrations == "vtk"} {
puts $fileId "2				! writing out concentrations in vtk using header:"
puts $fileId "\"$concentration_filename\""
for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {
puts $fileId "\"$concentration_header($jj)\""
}
} elseif {$write_concentrations == "no"} { puts $fileId "0				! not writing out concentrations"}

if {$write_out_particles == "yes"} {puts $fileId "1				! writing out particle locations to file:"
} else { puts $fileId "0				! not writing out particle locations"}
puts $fileId "\"$particle_out_file\""

#AMD - add for final output of particles
if {$write_out_end_particles == "yes"} {puts $fileId "1                        ! writing out final particle locations to file:"
} else {puts $fileId "0                                ! not writing out final particle locations"}
puts $fileId "\"$particle_out_end_file\""

if {$write_out_moments == "yes"} {puts $fileId "1				! writing out moments to file:"
} else { puts $fileId "0				! not writing out moments"}
puts $fileId "\"$moment_out_file\""

if {$simulation_mode == "forward"} {puts $fileId "0				! forward simulation"
} else { puts $fileId "1				! backward simulation"}

puts $fileId "$diffusivity						! molucular diffusivity,  \[m**2/day\]"

puts $fileId "$npmax						! maximum number of particles"
puts $fileId "$give_up						! maximum number of particle steps per time loop"
puts $fileId "$epsilon						! epsilon/zero value drop tolerance"


puts $fileId "$number_wells						! number of wells \[x,y,zt,zb,q\]"
for {set ii 1} {$ii <= $number_wells} {incr ii 1} {
puts $fileId "$well_x_location($ii) , $well_y_location($ii) , $well_screen_top($ii) , $well_screen_bottom($ii) , $well_pumping_rate($ii) "
}

if {$immobile_pause == "no"} {puts $fileId "0                   !running traditional SLIM"
} else { puts $fileId "1                                 !running time fractional SLIM"
puts $fileId "$lambda                                           !lambda"
puts $fileId "$gamma                                            !mass transfer rate coefficient (gamma)"
puts $fileId "$beta                                           !immobile/mobile mass ratio"}

for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {
puts $fileId "\"$plane_out_file($jj)\""
}


puts $fileId "$number_planes						! number of planes \[x/y/z,coord\]"

for {set ii 1} {$ii <= $number_planes} {incr ii 1} {
if {$plane_xyz($ii) == "x"} {puts $fileId "1				! plane in x"}
if {$plane_xyz($ii) == "y"} {puts $fileId "2				! plane in y"}
if {$plane_xyz($ii) == "z"} {puts $fileId "3				! plane in z"}
puts $fileId "$plane_location($ii)					! location of plane \[m\]"
}

for {set jj 1} {$jj <= $num_constituents} {incr jj 1} {

if {$slim_initial_condition_type($jj) == "well"} {
puts $fileId "1							! Well-Type IC"
puts $fileId "$well_ic_x($jj)				! Well IC: X \[m\] "
puts $fileId "$well_ic_y($jj)				! 	Y \[m\]"
puts $fileId "$well_ic_bot($jj)				! Zlower of well screen \[m\]"
puts $fileId "$well_ic_top($jj)				! Zupper of well screen \[m\]"
puts $fileId "$well_ic_pump_rate($jj)				! pumping rate \[vol/day\] "
puts $fileId "$well_ic_num_part_per($jj)				! number of particles per node "
puts $fileId "$well_ic_mass_per_part($jj)                        !mass per particle"
}

if {$slim_initial_condition_type($jj) == "stream"} {
puts $fileId "6							! StreamLine IC"
puts $fileId "$well_ic_x($jj)				! Well IC: X \[m\] "
puts $fileId "$well_ic_y($jj)				! 	Y \[m\]"
puts $fileId "$well_ic_bot($jj)				! Zlower of well screen \[m\]"
puts $fileId "$well_ic_top($jj)				! Zupper of well screen \[m\]"
puts $fileId "$well_ic_pump_rate($jj)				! pumping rate \[vol/day\] "
puts $fileId "$well_ic_num_part_per($jj)				! number of particles per node "
}

if {$slim_initial_condition_type($jj) == "pulse"} {
puts $fileId "2							! Pulse-Type IC"
puts $fileId "$slim_pulseIC_Xlower($jj)				! Pulse IC: Xlower \[m\] "
puts $fileId "$slim_pulseIC_Xupper($jj)                              ! Xupper \[m\] "
puts $fileId "$slim_pulseIC_Ylower($jj)				! Ylower \[m\]"
puts $fileId "$slim_pulseIC_Yupper($jj)                             ! Yupper \[m\]"
puts $fileId "$slim_pulseIC_Zlower($jj)				! Zlower \[m\]"
puts $fileId "$slim_pulseIC_Zupper($jj)				! Zupper \[m\]"
puts $fileId "$slim_pulseIC_Initial_Concentration($jj)			! Co \[mass/vol**3\]"
puts $fileId "$slim_pulseIC_number_particles($jj)				! num particles total "
puts $fileId "$slim_pulseIC_decay_rate($jj)				! pulse decay rate "
puts $fileId "$slim_pulseIC_decay_time($jj)				! pulse decay time \[days\] "
puts $fileId "$slim_pulseIC_decay_timesteps($jj)				! decay timesteps \[time/timesteps\] "
}
if {$slim_initial_condition_type($jj) == "ind_file"} {
puts $fileId "3							! Indicator file-Type IC"
puts $fileId "\"$slim_ic_ind_file($jj)\""
puts $fileId "$slim_ic_ind_cat($jj)						! indicator cat"
puts $fileId "$slim_Initial_Concentration($jj)			! Co \[mass/vol**3\]"
puts $fileId "$slim_number_particles($jj)				! num particles total "
}

if {$slim_initial_condition_type($jj) == "ind_mult"} {
puts $fileId "4							! Mulitple Indicator; reading in file-Type IC"
puts $fileId "\"$slim_ic_ind_file($jj)\""
}

if {$slim_initial_condition_type($jj) == "none"} {
puts $fileId "0							! No IC"
}

if {$slim_initial_condition_type($jj) == "cont"} {
puts $fileId "5                                                 ! Continued simulation"
puts $fileId "\"$slim_contIC_file\"                         ! IC input file"
puts $fileId "$slim_contIC_starttime                    ! Starttime"
puts $fileId "$slim_contIC_startvalue                   ! Startvalue"
}
}

# particle splitting?
if {$part_split == "yes"} {puts $fileId "1							! Particle Splitting"
puts $fileId "$min_conc						! min conc"
} else { puts $fileId "0							! Fixed Number of Particles" }
# temporal averaging?
if {$temp_averaging == "yes"} {puts $fileId "1							! Temporal Averaging"
} else { puts $fileId "0							! Not temporally avging for concentration" }
puts $fileId "$vel_nskip						! number of timesteps to skip"
puts $fileId " "
puts $fileId " "

puts $fileId " "
puts $fileId " "

close $fileId
#
#  Run SLIM 
#
set fileId [open "slim_dummy_file.par" w 0600] 
puts $fileId "$parameter_file"
close $fileId

exec $slim_executable  < slim_dummy_file.par > slim.out.txt

