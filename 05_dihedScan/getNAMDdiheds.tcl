
# Purpose: measure actual angle from dihedral constraint
#   in order to subtract the restraint energy from final energy.
# Usage: vmdt -e file.tcl

set hdir /work/cluster/limvt/hv1/05_dihedScan/2_tautomer/angles-md
cd $hdir
set dirlist [glob */ ]

set f [open diheds-from-coor.dat w] 

set ind1 1
set ind2 5
set ind3 6
set ind4 7


set tmpmolid 0
foreach i $dirlist {

    # get file names
    set dir [string trim $i "/"]
    set coorF [glob ${dir}/*coor]

    # load files into vmd
    mol new ${hdir}/../gbi-tautomer2.psf
    mol addfile ${hdir}/${dir}/gbi2-${dir}.coor

    # measure and write to file
    set angle [measure dihed [list [list $ind1 $tmpmolid] [list $ind2 $tmpmolid] [list $ind3 $tmpmolid] [list $ind4 $tmpmolid] ]]
    puts $f "$dir $angle"

    mol delete all
    set tmpmolid [expr $tmpmolid +1 ]
}
close $f
exit

