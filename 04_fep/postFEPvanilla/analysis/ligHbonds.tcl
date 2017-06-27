
# usage `vmd -e ligHbonds.tcl -args t1_d112e`
# not vmdt because plugins don't get loaded

set key [lindex $argv 0]

mol new $key.psf
mol addfile $key.dcd first 0 last -1 step 100 waitfor all 

hbonds -sel1 [atomselect top protein] -sel2 [atomselect top "resname GBI1"] -writefile yes -upsel yes -frames all -dist 3.5 -ang 40 -plot yes -log ${key}_hbonds.log -writefile yes -outfile ${key}_hbonds.dat -polar yes -DA both -type unique -detailout ${key}_hbonds-details.dat
