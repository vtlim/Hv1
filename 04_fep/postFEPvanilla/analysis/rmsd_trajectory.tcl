#RMSD script to run in VMD

# vmdt -e file.tcl -args pdbref withGBI inpsf indcd
#
set pdbref [lindex $argv 0]
set withGBI [lindex $argv 1]

set inpsf [lindex $argv 2]
set indcd [lindex $argv 3]

# open the trajectory to run the RMSD calculation.
mol new $pdbref
set molID 1

mol new $inpsf
mol addfile $indcd waitfor all $molID

# file to output data for plotting
set outDataFile [open rmsd_byChain.dat w]
if {$withGBI == 1 || $withGBI == 2} {
  puts $outDataFile "# Win | helix backbone rmsd | TM helix 1 | TM helix 2 | TM helix 3 | TM helix 4 | 2GBI rmsd (Angstroms)"
} else {
  puts $outDataFile "# Win | helix backbone rmsd | TM helix 1 | TM helix 2 | TM helix 3 | TM helix 4 (Angstroms)"
}

# set frame 0 as the reference
set refprot [atomselect 0 "protein and backbone and {{resid 100 to 125} or {resid 134 to 160} or {resid 168 to 191} or {resid 198 to 220}}" frame 0]
set refprot1 [atomselect 0 "protein and backbone and {{resid 100 to 125}}" frame 0]
set refprot2 [atomselect 0 "protein and backbone and {{resid 134 to 160}}" frame 0]
set refprot3 [atomselect 0 "protein and backbone and {{resid 168 to 191}}" frame 0]
set refprot4 [atomselect 0 "protein and backbone and {{resid 198 to 220}}" frame 0]
if {$withGBI == 1} {set refgbi [atomselect 0 "segname GBI1 and noh" frame 0] }
if {$withGBI == 2} {set refgbi [atomselect 0 "segname GBI2 and noh" frame 0] }

# selections being compared.
set compprot [atomselect $molID "protein and backbone and {{resid 100 to 125} or {resid 134 to 160} or {resid 168 to 191} or {resid 198 to 220}}"]
set compprot1 [atomselect $molID "protein and backbone and {{resid 100 to 125}}"]
set compprot2 [atomselect $molID "protein and backbone and {{resid 134 to 160}}"]
set compprot3 [atomselect $molID "protein and backbone and {{resid 168 to 191}}"]
set compprot4 [atomselect $molID "protein and backbone and {{resid 198 to 220}}"]
if {$withGBI == 1} {set compgbi [atomselect $molID "segname GBI1 and noh"]}
if {$withGBI == 2} {set compgbi [atomselect $molID "segname GBI2 and noh"]}

# calc rmsd for each frame.
set num_steps [molinfo $molID get numframes]
for {set frame 0} {$frame < $num_steps} {incr frame} {
    # get frame of interest
    $compprot frame $frame
    $compprot1 frame $frame
    $compprot2 frame $frame
    $compprot3 frame $frame
    $compprot4 frame $frame
    if {$withGBI == 1 || $withGBI == 2} {$compgbi frame $frame}

    # compute transformation & do alignment
    set trans_mat [measure fit $compprot $refprot]
    $compprot move $trans_mat
    if {$withGBI == 1 || $withGBI == 2} {$compgbi move $trans_mat}

    # compute the RMSD
    set rmsdprot [measure rmsd $compprot $refprot]
    set rmsdprot1 [measure rmsd $compprot1 $refprot1]
    set rmsdprot2 [measure rmsd $compprot2 $refprot2]
    set rmsdprot3 [measure rmsd $compprot3 $refprot3]
    set rmsdprot4 [measure rmsd $compprot4 $refprot4]
    if {$withGBI == 1 || $withGBI == 2} {
      set rmsdgbi [measure rmsd $compgbi $refgbi]
      lappend listgbi $rmsdgbi
    }

    # write info
    if {$withGBI == 1 || $withGBI == 2} {
      puts $outDataFile "$frame \t $rmsdprot \t $rmsdprot1 \t $rmsdprot2 \t $rmsdprot3 \t $rmsdprot4 \t $rmsdgbi"
    } else {
      puts $outDataFile "$frame \t $rmsdprot \t $rmsdprot1 \t $rmsdprot2 \t $rmsdprot3 \t $rmsdprot4"
    }
}
exit
