
# Purpose: Loop over all residues in the protein and print out non-neutral ones.
# Usage: Load in PSF and PDB, then in vmd type "source chargedRes.tcl"
#   * /beegfs/DATA/mobley/limvt/hv1/04_fep/sandbox_chgScript1/chargedRes.tcl
# Author: Victoria Lim
# Note: End terminal residues are considered with capping domain
#   * e.g., N terminal with PHE will look like +1 bc of the NH3+
#   * e.g., C terminal with ARG will look like 0 bc +1 of ARG cancels with 0 of COO-


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

set sel [atomselect top protein]
set chg [vecsum [$sel get charge]]
# alternatively: measure sumweights $sel weight charge
puts "\nTotal charge on protein is: $chg"
