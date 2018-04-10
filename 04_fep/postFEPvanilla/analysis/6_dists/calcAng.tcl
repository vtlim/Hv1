
# Purpose: Measure dihedral angle of atom1-atom2-atom3-atom4.
# Usage: vmd -dispdev none -e calcDist.tcl -args input.psf input.dcd atom-selection1 atom-selection2 atom-selection3 atom-selection4
# Note: To specify selection criteria, separate words with commas, such as: protein,and,resid,112
#
# ========================== Variables ========================= #

set inpsf [lindex $argv 0]
set indcd [lindex $argv 1]
set pre0 [lindex $argv 2]
set pre1 [lindex $argv 3]
set pre2 [lindex $argv 4]
set pre3 [lindex $argv 5]

# processing of selection arguments
set sel0 [split $pre0 {,}]
set sel1 [split $pre1 {,}]
set sel2 [split $pre2 {,}]
set sel3 [split $pre3 {,}]

set skip 10   ;# use 1 for no skip
set step0 0

# =============================================================== #

lappend auto_path /data12/cmf/limvt/tempotools/libs
lappend auto_path /home/limvt/Documents/tempotools/libs
package require tempoUserVMD
package require pbctools

# =============================================================== #


# define output file
set outDataFile [open measAng.dat w]
puts $outDataFile "# Input PSF: $inpsf\n# Input DCD, skip $skip: $indcd\n "
puts $outDataFile "# Atom 1: $sel0\n# Atom 2: $sel1\n# Atom 3: $sel2\n# Atom 3: $sel3"
puts $outDataFile "\n# Frame | Dihedral angle (degrees)"

# load files
mol new $inpsf
mol addfile $indcd type {dcd} first $step0 last -1 step $skip waitfor -1
pbc wrap -compound fragment -center com -centersel "protein and noh" -all 

# set atom selections
set atom1 [[atomselect top "$sel0"] get index] 
set atom2 [[atomselect top "$sel1"] get index]
set atom3 [[atomselect top "$sel2"] get index]
set atom4 [[atomselect top "$sel3"] get index]

# take measurements
set dihedlist [measure dihed [list $atom1 $atom2 $atom3 $atom4] frame all]

# write output
for {set i 0} {$i < [llength $dihedlist]} {incr i} {
    set ang [lindex $dihedlist $i]
    if {$ang < 0} {set ang [expr {360+$ang}]}
    puts $outDataFile "$i\t$ang" 
}
close $outDataFile

exit

