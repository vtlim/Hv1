### Tcl script to perform contact analyis on Hv1 + GBI trajectory
### Measure distances between GBI moities and important SF residues
### (See Hong et al. PNAS 2014)
### Usage: vmd -dispdev none -e file.tcl
### Author: Andrew Geragotelis, Tobias Lab @ UCI

# load trajectory
mol new /tw/ageragot/hv1/mdruns-hv1+gbi/hv1-config-15183_gbi-pose-04/hHv1_open_wGBI.psf
#mol addfile /tw/ageragot/hv1/mdruns-hv1+gbi/hv1-config-15183_gbi-pose-04/npt01.dcd waitfor all
mol addfile /pub/limvt/hv1/02_configs/15183_04/npt02.dcd waitfor all
mol addfile /pub/limvt/hv1/02_configs/15183_04/npt03.dcd waitfor all
mol addfile /pub/limvt/hv1/02_configs/15183_04/npt04.dcd waitfor all
mol addfile /pub/limvt/hv1/02_configs/15183_04/npt05.dcd waitfor all
mol addfile /pub/limvt/hv1/02_configs/15183_04/npt06.dcd waitfor all

# open data output file for making plots
set outFile [open "hv1+gbi_contacts.dat" w]
puts $outFile "#Frame | F150-benzo dist | R211-guan dist | D112-imidazole dist | S181-imidazole dist | R211-imidazole dist (Angstroms) "

# get number of trajectory frames
set nframes [molinfo top get numframes]
set all [atomselect top all]

# set selections
# CZ on R211 and guan C on 2GBI.
set index1 [[atomselect top "protein and resid 211 and name CZ"] get index]
set index2 [[atomselect top "resname GBI1 and name C"] get index]
# D112 sidechain O and imidizole N on 2GBI.
set index3a [[atomselect top "protein and resid 112 and name OD1"] get index]
set index3b [[atomselect top "protein and resid 112 and name OD2"] get index]
set index4 [[atomselect top "resname GBI1 and name N3"] get index]
# S181 sidechain O and imidizole N on 2GBI.
set index5 [[atomselect top "protein and resid 181 and name OG"] get index] 
set index6 [[atomselect top "resname GBI1 and name N4"] get index]

# loop over dcd frames
for {set i 0} {$i < $nframes} {incr i} {
    $all frame $i

    # measure distance between F150 phenyl ring to GBI benzo ring.
    # use geometric center of benz rings
    set vec1 [measure center [atomselect top "protein and resid 150 and name CG CD2 CE2 CZ CE1 CD1" frame $i]]
    set vec2 [measure center [atomselect top "resname GBI1 and name C2 C3 C4 C5 C6 C7" frame $i]]
    set Benzdist [veclength [vecsub $vec1 $vec2]]

    # measure distance between R211 and GBI guanidine moity
    # use guan carbon atom for both selections
    set guandist [measure bond "$index1 $index2" frame $i]


    # measure distance between D112,R211,S181 and GBI imidazole ring
    # use closest =O on D112 and N on 2GBI.
    set aspimidDist [expr min([measure bond "$index3a $index4" frame $i], [measure bond "$index3b $index4" frame $i])]
    set serimidDist [measure bond "$index5 $index6" frame $i]
    set argimidDist [expr min([measure bond "$index1 $index4" frame $i], [measure bond "$index1 $index6" frame $i])]

    # print data to output file
    puts $outFile "$i  $Benzdist  $guandist  $aspimidDist  $serimidDist  $argimidDist "
    unset vec1 vec2 Benzdist guandist aspimidDist serimidDist argimidDist
}

close $outFile
exit
