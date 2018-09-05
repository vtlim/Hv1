

# ======================================================================
# Purpose: Create wrapped summary trajectory.
#
# Usage: vmdt -e sumtraj.tcl -args psffile outputName skip dcdfile1 lastframe1 [dcdfile2 lastframe2] [dcdfile3 lastframe3]
# Ex:    vmdt -e sumtraj.tcl -args in.psf  wrap10.dcd 10   npt01.dcd 5000      npt02.dcd 5000
#
# Last frame can be obtained by subtracting one from catdcd output
# Version: Oct 4 2017
# ======================================================================

# only write out specific selection; else comment this line
#set extract "protein or resname GBI2"
#set extract "protein"

# for dopbc
lappend auto_path /data12/cmf/limvt/tempotools/libs
package require tempoUserVMD

# read in command line arguments
set arglen [llength $argv]
set inpsf  [lindex $argv 0]
set outputName [lindex $argv 1]
set skip  [lindex $argv 2]

# ======================================================================

# load protein structure file
mol new $inpsf

set index 3     ;# counter for file inputs (0=psf,1=output,2=skip)
set findex 4    ;# counter for frame number as index+1 (TODO redundant)
while {$index < $arglen} {

    # get command line arguments from index
    set dcdfile [lindex $argv $index]
    set lastf [lindex $argv $findex]

    # load and wrap trajectory
    dopbc -file $dcdfile -frames 0:$skip:$lastf -ref protein
    incr index 2
    incr findex 2
}


# write output
if {! [info exists extract] } {
    animate write dcd $outputName waitfor all
    exit
} else {
    # get psf basename
    set fbasename [file rootname [file tail $outputName]]
    # write dcd. if want subset of frames, use something like:   beg 0 end 9
    animate write dcd $outputName       sel [atomselect top "$extract"] waitfor all
    animate write psf "$fbasename.psf"  sel [atomselect top "$extract"]
    animate write pdb "$fbasename.pdb"  sel [atomselect top "$extract"] beg 0 end 0
    exit
}




