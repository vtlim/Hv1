
# Purpose:  Open all lambda windows in FEP simulation for viewing/analyzing.
#
# Usage:    vmd -e file.tcl -args homedir withGBI
#           - homedir: full path ending with pose and mutation, e.g. /path/to/19415-10/a1_F150A
#           - withGBI: 0 noGBI, 1 taut1, 2 taut2
#
# Notes: 
#  - If not 250 frames (5 ns), edit dopbc line down below.
#  - Can change skip value of trajectory reading in the dopbc line too.
#  - For opening just forward (F), or just reverse (R), or both, check outer foreach loop.
#
# File locations
#  - /beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/view/openAllWindows.tcl
#  - /home/limvt/connect/greenplanet/goto-beegfs/hv1/04_fep/analysis/view/openAllWindows.tcl 
#


# ========================== Variables ========================= #

set homedir [lindex $argv 0]
set withGBI [lindex $argv 1]

# Text processing
set pose [lindex [split $homedir /] end-1]
set mut [lindex [split $homedir /] end]

# Encoded variables, subject to change
set num_wins 40   ;# number of windows for each of fwd and rev

# Find input psf, original coordinates.
set psf [glob -dir ${homedir}/00_main *psf]
set pdb [glob -dir ${homedir}/00_main *pdb]


# =============================================================== #

lappend auto_path /data12/cmf/limvt/tempotools/libs
lappend auto_path /home/limvt/Documents/tempotools/libs
package require tempoUserVMD

# load in psf and reference coords
mol addfile $pdb waitfor all
mol new $psf
set molID 1

# loop over directions, fwd and rev
foreach way [list "R"] {

    # get list of trajectories, add zero for nums<10
    set dcdlist {}
    for {set i 1} {$i <= $num_wins} {incr i} {
        if {$i < 10} {
          set win 0$i
        } else {
          set win $i
        }
        lappend dcdlist ${homedir}/FEP_${way}/lambda_${win}/alchemy${win}.dcd
    }

    # loop over trajectories
    foreach dcd $dcdlist {
        dopbc -file ${dcd} -frames 0:10:249 -mol $molID -ref protein ;# get last frame (assuming 250 total)
    }

}
 
