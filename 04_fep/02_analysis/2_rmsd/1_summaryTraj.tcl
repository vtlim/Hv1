# Purpose of the script: to take and wrap every nth frame for summary trajectory.
# Usage: change files, then "vmdt -e file.tcl"
# vim to replace all foo with bar: :%s/foo/bar/g
# Tips: 
#  - To get endFrames, first need to know number of frames [catdcd -num file.dcd]
#    Then specify frames like [-frames 499:1:500]


# ========================== Variables ========================= #

set skip 1

set pose "17041_19"
set mut "F150A"
set way "F"
set num_wins 40

set dir /data12/cmf/limvt/hv1/04_fep/${pose}/${mut}

# Edit file names here.
set psf ${dir}/00_main/17041_19-F150A.psf
set out ${dir}/02_analysis/2_rmsd/finalFrames-${way}.dcd

# Don't forget to edit which frame selections you want (dopbc line)!


# =============================================================== #

lappend auto_path /data12/cmf/limvt/tempotools/libs
package require tempoUserVMD

set dcdlist {}

mol new $psf

if {[file exists $out]} {
    puts "\n\nError, output file already exists: $out\n\n"
    exit
}


for {set i 1} {$i <= $num_wins} {incr i} {
    # process numbers < 10 to add a zero
    if {$i < 10} {
        set win 0$i
    } else {
        set win $i
    }

    # get list of trajectories
    lappend dcdlist ${dir}/FEP_${way}/lambda_${win}/alchemy${win}.dcd
}

# add the first frame of the first dcd as reference
dopbc -file [lindex $dcdlist 0] -frames 0:1:0 -ref protein


foreach dcd $dcdlist {
    #dopbc -file ${dcd} -frames 0:$skip:-1 -ref protein ;# get every skipth frame
    dopbc -file ${dcd} -frames 249:1:250 -ref protein ;# get last? frame
    #dopbc -file ${dcd} -frames 0:1:1 -ref protein ;# get first frame of each window
}

animate write dcd $out waitfor all

set fo [open "${dir}/02_analysis/2_rmsd/README" "a"] 
puts $fo "\n\nDetails for summary trajectory:"
puts $fo " - Name:\t$out"
puts $fo " - Date:\t[clock format [clock seconds]]"
puts $fo " - Skip:\t$skip"
puts $fo " - System:\t$pose/$mut/FEP_$way"
puts $fo " - PSF file:\t$psf"
close $fo
exit

