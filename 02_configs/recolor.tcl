# ligand can be colored by name, by type, or by element
# this script sets the element name so I can make use of all three
# ex. name for side chains, type for t2, and elem for t1

set sel1 [atomselect top "resname GBI1 and hydrogen"]
$sel1 set element H

set sel2 [atomselect top "resname GBI1 and nitrogen"]
$sel2 set element N

set sel3 [atomselect top "resname GBI1 and carbon"]
$sel3 set element C

