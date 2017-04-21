
# Purpose: 

# Usage: vmd -e file.tcl -args homedir withGBI
#   - homedir: full path ending with pose and mutation, e.g. /path/to/19415-10/a1_F150A
#   - withGBI: 0 noGBI, 1 taut1, 2 taut2
#

# ========================== Variables ========================= #

set homedir [lindex $argv 0]
set withGBI [lindex $argv 1]



# Text processing
set pose [lindex [split $homedir /] end-1]
set mut [lindex [split $homedir /] end]

# Encoded variables, subject to change
# If not 250 frames (5 ns), edit dopbc line down below.
set num_wins 40   ;# number of windows for each of fwd and rev



# Find input psf, original coordinates.
set psf [glob -dir ${homedir}/00_main *psf]
set pdb [glob -dir ${homedir}/00_main *pdb]


# Write out summary trajectory
# NOT YET IMPLEMENTED,
#  - take in skip value for reading in dcd
#  - don't delete molID 
set writeSumTraj 0




# =============================================================== #

lappend auto_path /data12/cmf/limvt/tempotools/libs
lappend auto_path /home/limvt/Documents/tempotools/libs
package require tempoUserVMD

# load in psf and reference coords
mol addfile $pdb waitfor all
mol new $psf
set molID 1

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

    foreach dcd $dcdlist {
        dopbc -file ${dcd} -frames 0:1:249 -mol $molID -ref protein ;# get last frame (assuming 250 total)
    }





}
 
