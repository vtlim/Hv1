

# Usage: vmdt -e trajComponent.tcl -args file.psf file.dcd

# Calculate some component of the energy for selected atoms in trajectory.

# ====================================================================================


# open data output file for making plots
set outFile [open "vdw-2gbi-protlip.dat" w]
puts $outFile "# Frame\tV_vdw "


set psf [lindex $argv 0]
set dcd [lindex $argv 1]

mol new $psf
mol addfile $dcd first 0 last -1  waitfor all

set lig [[atomselect top "resname GBI1"] get index] 

set env [[atomselect top "(lipid and resid 149) or (protein and within 6 of resname GBI1:"] get index]  ;# imidazolic H
set all [atomselect top all]

set nframes [molinfo 0 get numframes]
for {set i 0} {$i < $nframes} {incr i} {
    $all frame $i

    # measure the dihedral angle
    # not yet working... set dene [measure energy dihed "$a1 $a2 $a3 $a4" frame $i]
    set eene [measure energy elect "$h2 $h3" q1 0.29 q2 0.42 frame $i]

    # print data to output file
    puts $outFile "$i\t$eene"
    unset eene
}

close $outFile
exit
