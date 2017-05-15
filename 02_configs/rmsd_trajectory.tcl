
# Purpose: Measure backbone RMSD along trajectory using VMD.
# Usage: vmd -dispdev none -e file.tcl -args inpsf indcd withGBI
#   - withGBI: 0 noGBI, 1 taut1, 2 taut2
#
# Authors: Andrew Geragotelis, Victoria Lim | Tobias Lab @ UCI
# TODO: modify to take in a list of dcd, split by commas, and add by mol addfile



# ========================== Variables ========================= #

set inpsf [lindex $argv 0]
set indcd [lindex $argv 1]
set withGBI [lindex $argv 2]



# specify ligand
if {$withGBI == 1} {
    set GBI "resname GBI1 and noh"
} elseif {$withGBI == 2} {
    set GBI "resname GBI2 and noh"
} else {
    puts "Specify a valid version of 2GBI ligand."
    exit
}

# open the trajectory to run the RMSD calculation.
mol new $inpsf
#mol addfile $indcd type {dcd} first 0 last -1 step 1 waitfor -1
mol addfile $indcd type {dcd} first 0 last -1 step 20 waitfor -1

# file to output data for plotting
set outDataFile [open rmsd.dat w]
puts $outDataFile "#Frame | helix backbone rmsd | Res112,150,181,211 rmsd | 2GBI rmsd (Angstroms)"

# set frame 0 as the reference
set refprot [atomselect top "protein and backbone and {{resid 98 to 126} or {resid 134 to 160} or {resid 169 to 191} or {resid 198 to 221}}" frame 0]
set refres [atomselect top "protein and resid 150 112 211 181" frame 0]
set refgbi [atomselect top "$GBI" frame 0]

# selections being compared.
set compprot [atomselect top "protein and backbone and {{resid 98 to 126} or {resid 134 to 160} or {resid 169 to 191} or {resid 198 to 221}}"]
set compgbi [atomselect top "$GBI"]
set compres [atomselect top "protein and resid 150 112 211 181"]
set all [atomselect top all]

# calc rmsd for each frame.
set num_steps [molinfo 0 get numframes]
for {set frame 0} {$frame < $num_steps} {incr frame} {

    # get frame of interest
    $all frame $frame
    $compprot frame $frame
    $compres frame $frame
    $compgbi frame $frame

    # compute transformation
    set trans_mat [measure fit $compprot $refprot]

    # do the alignment. move all else sidechains won't get moved w/refprot
    $all move $trans_mat

    # compute the RMSD
    set rmsdprot [measure rmsd $compprot $refprot]
    set rmsdres [measure rmsd $compres $refres]
    set rmsdgbi [measure rmsd $compgbi $refgbi]

    # print to file
    puts $outDataFile "$frame \t $rmsdprot \t $rmsdres \t $rmsdgbi"
}

close $outDataFile
exit
