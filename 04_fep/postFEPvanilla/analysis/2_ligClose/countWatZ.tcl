
# Purpose: Count number of waters from -z0 to +z1 coordinate.

# Usage: vmdt -e file.tcl -args inpsf indcd selection_z0 selection_z1
#   - Selection for z0 should be the lower z plane
#   - Selection for z1 should be the upper z plane
#
# Modified from: http://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/23723.html
# Specify selection with no spaces, but use commas for multiple words: protein,and,resid,211


# ========================== Variables ========================= #

set inpsf [lindex $argv 0]
set indcd [lindex $argv 1]
set pre_z0 [lindex $argv 2]
set pre_z1 [lindex $argv 3]

set watlist [list]
set sel_z0 [split $pre_z0 {,}]
set sel_z1 [split $pre_z1 {,}]

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
set outDataFile [open waters-in-zrange.dat w]
puts $outDataFile "# Frame | number of waters bt Z crds of two sels"

# load files
mol new $inpsf
mol addfile $indcd type {dcd} first 0 last -1 step 10 waitfor -1
pbc wrap -compound fragment -center com -centersel "protein" -all 

# get reference Z coordinates
set z0 [[atomselect top $sel_z0] get {z}]
set z1 [[atomselect top $sel_z1] get {z}]
puts $outDataFile "# Selection 1 (z=$z0): $sel_z0\n# Selection 2 (z=$z1): $sel_z1"

# loop over frames
set frames [molinfo top get numframes]
set wats [atomselect top "noh and waters and (z<$z1 and z>$z0) and (x<15 and x>-20) and (y<20 and y>-10)"]
for {set i 0} {$i < $frames} {incr i} {
    $wats frame $i
    $wats update
    set num [$wats num]
    lappend watlist $num
    puts $outDataFile "$i $num"
}

$wats delete
puts $outDataFile "# --- Average over traj: [average $watlist]"
close $outDataFile


exit
