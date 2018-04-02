
# Purpose: Count number of waters within x Angstrom of specified selection.

# Usage: vmdt -e file.tcl -args inpsf indcd selection
#
# Modified from: http://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/23723.html
# Specify selection with no spaces, but use commas for multiple words: protein,and,resid,211


# ========================== Variables ========================= #

set inpsf [lindex $argv 0]
set indcd [lindex $argv 1]
set prelig [lindex $argv 2]

set dist 4
set watlist [list]
set lig [split $prelig {,}]

# =============================================================== #

proc average L {
    expr ([join $L +])/[llength $L].
}

lappend auto_path /data12/cmf/limvt/tempotools/libs
lappend auto_path /home/limvt/Documents/tempotools/libs
package require tempoUserVMD
package require pbctools

# =============================================================== #


# define output file
set outDataFile [open waters-within-${dist}.dat w]
puts $outDataFile "# Frame | number of waters (noh) within $dist Angstroms of $lig (noh)"

# load files
mol new $inpsf
mol addfile $indcd type {dcd} first 0 last -1 step 10 waitfor -1
#mol addfile $indcd type {dcd} first 0 last 100 step 1 waitfor -1

# specify ligand
set ligand "$lig and noh"
pbc wrap -compound fragment -center com -centersel "$ligand" -all 

# loop over frames
set frames [molinfo top get numframes]
set cw [atomselect top "water and oxygen within $dist of $ligand"]
for {set i 0} {$i < $frames} {incr i} {
    $cw frame $i
    $cw update
    set num [$cw num]
    lappend watlist $num
    puts $outDataFile "$i $num"
}

$cw delete
puts $outDataFile "# --- Average over traj: [average $watlist]"
close $outDataFile



exit

