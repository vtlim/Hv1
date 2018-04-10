

# Purpose: Process a NAMD .fepout file and load frames of the corresponding trajectory
#   when the work meets a certain threshold. Manually set the cutoff before using this script!
# Usage: vmdt -e file.tcl -args in.fepout in.psf in.dcd outputBase
# Example: vmdt -e extractWork.tcl -args ../../FEP_R/lambda_01/alchemy01.fepout ../../00_main/18629-19_R211S.psf ../../FEP_R/lambda_01/alchemy01.dcd R01 
# TODO: add functionality for multiple cutoffs (and turn code into functions)

# ================================================================================================

set work_cutoff 2.5     ; # cutoff value from which to separate trajectory
set jump_dcd 10000      ; # how many timesteps between one frame and the next in the DCD (DCDfreq)

# ================================================================================================


# set variables from arguments
set fepout [lindex $argv 0]
set mypsf [lindex $argv 1]
set mydcd [lindex $argv 2]
set outbase [lindex $argv 3]

# read in fepout file and convert lines to list. modified from:
# https://stackoverflow.com/questions/18219420/how-to-convert-a-file-to-a-tcl-list
set f [open $fepout "r"]
foreach line [split [read $f] \n] {
    # skip lines starting with #
    if {[string match "#*" $line]} {
        continue
    }
    lappend times [lindex $line 1]
    lappend works [lindex $line 6]
}
close $f

# remove null items from lists
set times [lsearch -all -inline -not -exact $times {}]
set works [lsearch -all -inline -not -exact $works {}]

# bin each timestep from work values based on work_cutoff
set bin1 {}
set bin2 {}
set num_snapshots [llength $works]
for {set i 0} {$i < $num_snapshots} {incr i} {
    if {[lindex $works $i] <= $work_cutoff} {
        lappend bin1 [lindex $times $i]
    } else {
        lappend bin2 [lindex $times $i]
    }
}

# write the binned values to file
set outfile [open "timesteps_work_$outbase.dat" w]
puts $outfile "# Analysis of NAMD .fepout file: $fepout"
puts $outfile "# Cutoff of dE values: $work_cutoff"
puts $outfile "# [llength $bin1] time steps below cutoff:\n$bin1"
puts $outfile "# [llength $bin2] time steps above cutoff:\n$bin2"
close $outfile

# extract/convert binned timesteps from .fepout to timesteps of dcd
set steps1 {}
set steps2 {}
foreach elem $bin1 {
    lappend steps1 [expr {$elem * 1.0 / $jump_dcd}]
}
foreach elem $bin2 {
    lappend steps2 [expr {$elem * 1.0 / $jump_dcd}]
}

# load in psf and dcd, by each set of binned values
# BIN 1
mol new $mypsf
foreach step $steps1 {
    if {$step == floor($step)} { ;# only consider whole number steps
        # convert to zero-based indices then read in DCD frame
        set i [expr { int($step - 1.0) }]
        puts $i
        animate read dcd $mydcd beg $i end $i waitfor all 0
    }
}
animate write dcd "traj_${outbase}_bin1.dcd" waitfor all 0

# BIN 2
mol new $mypsf
foreach step $steps2 {
    if {$step == floor($step)} { ;# only consider whole number steps
        # convert to zero-based indices then read in DCD frame
        set i [expr { int($step - 1.0) }]
        puts $i
        animate read dcd $mydcd beg $i end $i waitfor all 1
    }
}
animate write dcd "traj_${outbase}_bin2.dcd" waitfor all 1





exit
