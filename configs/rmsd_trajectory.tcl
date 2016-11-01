### RMSD script to run in VMD
### Usage: vmd -dispdev none -e file.tcl
### Author: Andrew Geragotelis, Tobias Lab @ UCI

# open the trajectory to run the RMSD calculation.
mol new /tw/ageragot/hv1/mdruns-hv1+gbi/hv1-config-15183_gbi-pose-04/hHv1_open_wGBI.psf
#mol addfile /tw/ageragot/hv1/mdruns-hv1+gbi/hv1-config-15183_gbi-pose-04/npt01.dcd waitfor all
mol addfile /pub/limvt/hv1/02_configs/15183_04/npt02.dcd waitfor all
mol addfile /pub/limvt/hv1/02_configs/15183_04/npt03.dcd waitfor all
mol addfile /pub/limvt/hv1/02_configs/15183_04/npt04.dcd waitfor all
mol addfile /pub/limvt/hv1/02_configs/15183_04/npt05.dcd waitfor all
mol addfile /pub/limvt/hv1/02_configs/15183_04/npt06.dcd waitfor all

# file to output data for plotting
set outDataFile [open rmsd_withStart.dat w]
puts $outDataFile "#Frame | helix backbone rmsd | Res112,150,181,211 rmsd | 2GBI rmsd (Angstroms)"

# set frame 0 as the reference
set refprot [atomselect top "protein and backbone and {{resid 98 to 126} or {resid 134 to 160} or {resid 169 to 191} or {resid 198 to 221}}" frame 0]
set refres [atomselect top "protein and resid 150 112 211 181" frame 0]
set refgbi [atomselect top "segname GBI1 and noh" frame 0]

# selections being compared.
set compprot [atomselect top "protein and backbone and {{resid 98 to 126} or {resid 134 to 160} or {resid 169 to 191} or {resid 198 to 221}}"]
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
    set rmsdres [measure rmsd $compres $refres]
    set rmsdgbi [measure rmsd $compgbi $refgbi]
    # print to file
    puts $outDataFile "$frame \t $rmsdprot \t $rmsdres \t $rmsdgbi"
}
exit
