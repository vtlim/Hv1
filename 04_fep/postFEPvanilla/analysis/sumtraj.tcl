

# Purpose: Create wrapped summary trajectory.
# Usage: vmdt -e sumtraj.tcl -args psffile skip outputName dcdfile1 lastframe1 [dcdfile2 lastframe2] [dcdfile3 lastframe3]
# Last frame can be obtained by subtracting one from catdcd output
# !! MULTIPLE DCD'S HAVE NOT BEEN TESTED AND VERIFIED !!

# for dopbc
package require tempoUserVMD

# read in command line arguments
set arglen [llength $argv]
set inpsf  [lindex $argv 0]
set skip  [lindex $argv 1]
set outputName [lindex $argv 2]


mol new $inpsf

set index 3
set findex 4
while {$index < $arglen} {
    set dcdfile [lindex $argv $index]
    set lastf [lindex $argv $findex]
    dopbc -file $dcdfile -frames 0:$skip:$lastf -ref protein
    incr index 2
    incr findex 2
}

animate write dcd $outputName waitfor all
