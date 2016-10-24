# Purpose of the script: to take and wrap every nth frame for summary trajectory.
# Usage: change files, then "vmdt -e file.tcl"
# vim to replace all foo with bar: :%s/foo/bar/g


# ========================== Variables ========================= #

set skip 1

set pose "15183_04"
set mut "F150A-noGBI"
set way "F"
set num_wins 40

set dir /data12/cmf/limvt/hv1/04_fep/${pose}/${mut}

# Edit file names at the very end.
set psf ${dir}/00_main/15183-F150A.psf
set out ${dir}/02_analysis/2_summaryTraj/endFrames.dcd


# =============================================================== #

lappend auto_path /data12/cmf/limvt/tempotools/libs
package require tempoUserVMD

set dcdlist {}

mol new $psf



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

foreach dcd $dcdlist {
    #dopbc -file ${dcd} -frames 0:$skip:-1 -ref protein
    dopbc -file ${dcd} -frames 249:1:250 -ref protein
}

animate write dcd $out waitfor all

set fo [open "${dir}/02_analysis/2_summaryTraj/README" "a"] 
puts $fo "\n\nDetails for summary trajectory:"
puts $fo " - Name:\t$out"
puts $fo " - Date:\t[clock format [clock seconds]]"
puts $fo " - Skip:\t$skip"
puts $fo " - System:\t$pose/$mut/FEP_$way"
puts $fo " - PSF file:\t$psf"
close $fo
exit

