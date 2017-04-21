
# Purpose: Get lambda 0 snapshots from fwd/rev, and lambda 1 snapshots from fwd/rev.

# Usage: vmdt -e file.tcl -args homedir withGBI
#   - homedir: full path ending with pose and mutation, e.g. /path/to/19415-10/a1_F150A
#   - withGBI: 0 noGBI, 1 taut1, 2 taut2
#
# Tips: 
#  - If not 250 frames (5 ns), edit dopbc line down below.
#  - To get endFrames, first need to know number of frames [catdcd -num file.dcd]
#    Then specify frames like [-frames 499:1:500]


# ========================== Variables ========================= #

set homedir [lindex $argv 0]
set withGBI [lindex $argv 1]

# Text processing
set pose [lindex [split $homedir /] end-1]
set mut [lindex [split $homedir /] end]

# Find input psf, original coordinates.
set psf [glob -dir ${homedir}/00_main *psf]
set pdb [glob -dir ${homedir}/00_main *pdb]
set fo [open "${homedir}/02_analysis/2_rmsd/README" "a"] 

# Archived.
#set skip 1


# =============================================================== #

lappend auto_path /data12/cmf/limvt/tempotools/libs
lappend auto_path /home/limvt/Documents/tempotools/libs
package require tempoUserVMD

# load in psf and reference coords
mol new $psf
mol addfile $pdb

# define output file
set outtraj ${homedir}/02_analysis/2_rmsd/endpoints.dcd   ;#

## check that summary trajectory doesn't already exist    
# if {[file exists $outtraj]} {
#     puts "\n\nError, output file already exists: $outtraj\n\n"
#     exit
#}

set dcdlist {}
lappend dcdlist ${homedir}/FEP_F/lambda_01/alchemy01.dcd ;# lambda0_fwd
lappend dcdlist ${homedir}/FEP_R/lambda_40/alchemy40.dcd ;# lambda0_rev
lappend dcdlist ${homedir}/FEP_F/lambda_40/alchemy40.dcd ;# lambda1_fwd
lappend dcdlist ${homedir}/FEP_R/lambda_01/alchemy01.dcd ;# lambda1_rev

set starts [list 50 249 249 50]
set ends [list 50 249 249 50]

# loop over specified lambda windows
foreach dcd $dcdlist sf $starts ef $ends {
    dopbc -file ${dcd} -frames $sf:1:$ef -ref protein ;# get last frame (assuming 250 total)
}

    
# Write out summary trajectory
animate write dcd $outtraj waitfor all


# Record details of summary trajectory in README 
puts $fo "\n\nDetails for summary trajectory:"
puts $fo " - Name:\t$outtraj"
puts $fo " - Date:\t[clock format [clock seconds]]"
#puts $fo " - Skip:\t$skip"
puts $fo " - System:\t$pose/$mut"
puts $fo " - PSF:\t\t$psf"
close $fo

exit
 
