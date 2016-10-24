# RMSD script to run in VMD
# Usage: vmdt -e file.tcl
# Adapted from Hv1 RMSD trajectory script from Andrew running poses

# open the trajectory to run the RMSD calculation.
mol new /data12/cmf/limvt/hv1/04_fep/15183_04/F150A-noGBI/00_main/15183-F150A.psf
mol addfile /data12/cmf/limvt/hv1/04_fep/15183_04/F150A-noGBI/02_analysis/2_summaryTraj/endFrames.dcd waitfor all
#mol addfile /pub/limvt/hv1/02_configs/15183_04/npt02.dcd waitfor all

# file to output data for plotting
set outDataFile [open rmsd_endFrames.dat w]
puts $outDataFile "#Frame | helix backbone rmsd | Res112,150,181,211 rmsd | 2GBI rmsd (Angstroms)"
#puts $outDataFile "#Frame | helix backbone rmsd | Res112,150,181,211 rmsd | 2GBI rmsd (Angstroms)"

# set frame 0 as the reference
set refprot [atomselect top "protein and backbone and {{resid 100 to 125} or {resid 134 to 160} or {resid 168 to 191} or {resid 198 to 220}}" frame 0]
set refres [atomselect top "protein and resid 150 112 211 181" frame 0]
set refgbi [atomselect top "segname GBI1 and noh" frame 0]

# selections being compared.
set compprot [atomselect top "protein and backbone and {{resid 100 to 125} or {resid 134 to 160} or {resid 168 to 191} or {resid 198 to 220}}"]
set compres [atomselect top "protein and resid 150 112 211 181"]
set compgbi [atomselect top "segname GBI1 and noh"]

set all [atomselect top all]

# calc rmsd for each frame.
set num_steps [molinfo 0 get numframes]
for {set frame 0} {$frame < $num_steps} {incr frame} {
    # get frame of interest
    $compprot frame $frame
    $compres frame $frame
    $compgbi frame $frame
    # compute transformation
    set trans_mat [measure fit $compprot $refprot]
    # do the alignment
    $compprot move $trans_mat
    $compres move $trans_mat
    $compgbi move $trans_mat
    # compute the RMSD
    set rmsdprot [measure rmsd $compprot $refprot]
#    set rmsdres [measure rmsd $compres $refres]
#    set rmsdgbi [measure rmsd $compgbi $refgbi]
    # print to file
    puts $outDataFile "$frame \t $rmsdprot"
#    puts $outDataFile "$frame \t $rmsdprot \t $rmsdres \t $rmsdgbi"
}
exit
