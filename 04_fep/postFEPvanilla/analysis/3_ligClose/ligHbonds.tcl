
# Purpose: Quantify number of hbonds from (protein or water) to (given selection).
# Usage:   vmdt -e ligHbonds.tcl -args selection,phrase,here inpsf mydcd1 mydcd2 [...]
# Example: vmdt -e ligHbonds.tcl -args protein,and,resid,211 file.psf traj1.dcd

package require pbctools
package require hbonds

# ========================== Variables ========================= #

set inskip 10

set prelig [lindex $argv 0]
set inpsf [lindex $argv 1]
for {set i 2} {$i < $argc} {incr i} {
    lappend dcdlist [lindex $argv $i]
}
set mylig [split $prelig {,}]

# =============================================================== #

# load system and wrap trajectory
mol new $inpsf
foreach dcd $dcdlist {    ;# maybe alter the first step to read in if FEP bc 50 frames equil
    mol addfile $dcd first 0 last -1 step $inskip waitfor all
}
pbc wrap -compound fragment -center com -centersel "$mylig" -all

# evaluate hbonds
hbonds -sel1 [atomselect top "protein or water"] -sel2 [atomselect top "$mylig"] -writefile yes -upsel yes -frames all -dist 3.5 -ang 40 -plot yes -log hbondsProtWat.log -writefile yes -outfile hbondsProtWat.dat -polar yes -DA both -type unique -detailout hbondsProtWat-details.dat

# append trajectory information to output
set outDataFile [open hbondsProtWat.log a]
puts $outDataFile "\n# Input PSF: $inpsf\n# Input DCD, skip $inskip: $dcdlist"
puts $outDataFile "# Hbonds selection: $mylig\n"
close $outDataFile

exit

