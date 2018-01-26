
# Purpose: Loop over all residues in the protein and print out non-neutral ones.
# Usage:
# Author: Victoria Lim


set sel [atomselect top protein]
set reslist [lsort -integer -unique [$sel get resid]]

foreach ri $reslist {
    set sel [atomselect top "protein and resid $ri"]
    set rn [lindex [[atomselect top "protein and resid $ri"] get resname] 0]
    set chg [vecsum [$sel get charge]]
    set rounded [expr round($chg)]
    if {$rounded != 0} {
        puts "$ri $rn\t$rounded"
    }

}

