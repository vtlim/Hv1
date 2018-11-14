
# Purpose: Prepare PDB file for REST2 simulations
# Usage: vmd -dispdev text -e prep_rest2_region.tcl -args input.pdb vmd,atom,selection output.pdb

# read in command line arguments
set inpdb  [lindex $argv 0]
set vmd_sel_pre [lindex $argv 1]
set outpdb  [lindex $argv 2]

# process vmd selection
set vmd_sel [split $vmd_sel_pre {,}]

# do the job
mol new $inpdb
set all [atomselect top all]
$all set occupancy 0.0
set sel [atomselect top "$vmd_sel"]
$sel set occupancy 1.0
animate write pdb $outpdb
exit

