
# Purpose: Measure distance of selection1 to selection[2,3]...
# Usage: vmd -dispdev none -e calcDist.tcl -args input.psf input.dcd atom-selection1 atom-selection2 [atom-selection3] [atom-selection4]
# Notes: 
#  - If calcDistances to match with colvars traj file, uncomment the step0 line below the line "set skip 10"
#  - If you get error: "measure center: bad weight sum, would cause divide by zero", double check your selection criteria.
#  - To specify selection criteria, separate words with commas, not spaces. ex: protein,and,resid,112
#
# ========================== Variables ========================= #

set inpsf [lindex $argv 0]
set indcd [lindex $argv 1]
set pre0 [lindex $argv 2]
set pre1 [lindex $argv 3]
set pre2 [lindex $argv 4]
set pre3 [lindex $argv 5]

# processing of other arguments
set sel0 [split $pre0 {,}]
set sel1 [split $pre1 {,}]
if { $argc >= 6 } {set sel2 [split $pre2 {,}]}
if { $argc >= 7 } {set sel3 [split $pre3 {,}]}

set skip 10                 ;# use 1 for no skip
#set step0 [expr $skip - 1] ;# to match start of colvars traj file 
set step0 0

# =============================================================== #

lappend auto_path /data12/cmf/limvt/tempotools/libs
lappend auto_path /home/limvt/Documents/tempotools/libs
package require tempoUserVMD
package require pbctools

# =============================================================== #


# define output file
set outDataFile [open measDist.dat w]
puts $outDataFile "# Input PSF: $inpsf\n# Input DCD, skip $skip: $indcd\n"
puts $outDataFile "# Selection 0: $sel0\n# Selection 1: $sel1"

if {[info exists sel2]} {
  puts $outDataFile "# Selection 2: $sel2"
}
if {[info exists sel3]} {
  puts $outDataFile "# Selection 3: $sel3"
}
puts $outDataFile "\n# Frame | Distance from sel0 to all other selections (Angstrom)"

# load files
mol new $inpsf
mol addfile $indcd type {dcd} first $step0 last -1 step $skip waitfor -1 ; # match start of colvars traj file
pbc wrap -compound fragment -center com -centersel "protein and noh" -all 

# set atom selections
set atom1 [[atomselect top "$sel0"] get index] 
set atom2 [[atomselect top "$sel1"] get index]
if {[info exists sel2]} {set atom3 [[atomselect top "$sel2"] get index]}
if {[info exists sel3]} {set atom4 [[atomselect top "$sel3"] get index]}

# take measurements
set dist1 [measure bond [list $atom1 $atom2] frame all]
if {[info exists sel2]} {set dist2 [measure bond [list $atom1 $atom3] frame all]}
if {[info exists sel3]} {set dist3 [measure bond [list $atom1 $atom4] frame all]}

# write output
for {set i 0} {$i < [llength $dist1]} {incr i} {
    if {[info exists dist3]} { puts $outDataFile "$i\t[lindex $dist1 $i]\t[lindex $dist2 $i]\t[lindex $dist3 $i]" } \
    elseif {[info exists dist2]} { puts $outDataFile "$i\t[lindex $dist1 $i]\t[lindex $dist2 $i]" } \
    else { puts $outDataFile "$i\t[lindex $dist1 $i]" }
}
close $outDataFile

exit

