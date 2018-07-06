
# usage `vmdt -e ligHbonds.tcl -args mypsf mydcd resname,GBI2`

package require pbctools
package require hbonds

set mypsf [lindex $argv 0]
set mydcd [lindex $argv 1]
set prelig [lindex $argv 2]

set mylig [split $prelig {,}]

mol new $mypsf
mol addfile $mydcd first 0 last -1 step 10 waitfor all 
pbc wrap -compound fragment -center com -centersel "$mylig" -all

hbonds -sel1 [atomselect top "protein or water"] -sel2 [atomselect top "$mylig"] -writefile yes -upsel yes -frames all -dist 3.5 -ang 40 -plot yes -log hbondsProtWat.log -writefile yes -outfile hbondsProtWat.dat -polar yes -DA both -type unique -detailout hbondsProtWat-details.dat

exit
