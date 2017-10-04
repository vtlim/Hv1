
# Purpose: Count number of waters within x Angstrom of ligand.

# Usage: vmdt -e file.tcl -args inpsf indcd withGBI
#   - withGBI: 0 noGBI, 1 taut1, 2 taut2
#
# Modified from: http://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/23723.html
#


# ========================== Variables ========================= #

set inpsf [lindex $argv 0]
set indcd [lindex $argv 1]
set withGBI [lindex $argv 2]

set dist 4
set watlist [list]
# =============================================================== #

proc average L {
    expr ([join $L +])/[llength $L].
}

lappend auto_path /data12/cmf/limvt/tempotools/libs
lappend auto_path /home/limvt/Documents/tempotools/libs
package require tempoUserVMD

# =============================================================== #


# define output file
set outDataFile [open waters-within-${dist}.dat w]
puts $outDataFile "# Frame | number of waters (noh) within $dist Angstrom of 2GBI (noh)"

# load files
mol new $inpsf
mol addfile $indcd type {dcd} first 0 last -1 step 1 waitfor -1
#mol addfile $indcd type {dcd} first 0 last 100 step 1 waitfor -1

# specify ligand
if {$withGBI == 1} {
    set ligand "resname GBI1 and noh"
} elseif {$withGBI == 2} {
    set ligand "resname GBI2 and noh"
} else {
    puts "Specify a valid version of 2GBI ligand."
    exit
}

# loop over frames
set frames [molinfo top get numframes]
for {set i 0} {$i < $frames} {incr i} {
    set cw [atomselect top "water and oxygen within $dist of $ligand"]
    $cw frame $i
    $cw update
    set num [$cw num]
    lappend watlist $num
    puts $outDataFile "$i $num"
    $cw delete
}

puts $outDataFile "# --- Average over traj: [average $watlist]"
close $outDataFile



exit

