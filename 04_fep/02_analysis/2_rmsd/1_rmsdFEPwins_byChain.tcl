
# Purpose: measure backbone RMSD of every last (or otherwise specified) frame 
#          of all fwd and rev FEP windows compared to initial seed coordinates (PDB).

# Usage: vmdt -e file.tcl -args homedir withGBI
#   - homedir: full path ending with pose and mutation, e.g. /path/to/19415-10/a1_F150A
#   - withGBI: 0 noGBI, 1 taut1, 2 taut2
#
# vim to replace all foo with bar: :%s/foo/bar/g
#
# Tips: 
#  - To get endFrames, first need to know number of frames [catdcd -num file.dcd]
#    Then specify frames like [-frames 499:1:500]
#
# Notes to self:
#  - don't do trans_mat on residues bc then they're doubly moved
#  - mol addfile yields same results as dopbc


# ========================== Variables ========================= #

set homedir [lindex $argv 0]
set withGBI [lindex $argv 1]



# Text processing
set pose [lindex [split $homedir /] end-1]
set mut [lindex [split $homedir /] end]

# Encoded variables, subject to change
# If not 250 frames (5 ns), edit dopbc line down below.
set num_wins 40   ;# number of windows for each of fwd and rev



# Find input psf, original coordinates.
set psf [glob -dir ${homedir}/00_main *psf]
set pdb [glob -dir ${homedir}/00_main *pdb]


# Write out summary trajectory
# NOT YET IMPLEMENTED,
#  - take in skip value for reading in dcd
#  - don't delete molID 
set writeSumTraj 0


# ================= User-defined functions ===================== #

proc average L {
    expr ([join $L +])/[llength $L].
}

proc rmsdTraj {molID withGBI} {
    # -------------------------------------------------------------
    # Parameters: VMD molID
    # Returns: the average RMSD of this trajectory for Hv1, 2GBI
    #          compared to initial coordinates
    # -------------------------------------------------------------


    # set selections from reference frame.
    set refprot [atomselect 0 "protein and backbone and {{resid 100 to 125} or {resid 134 to 160} or {resid 168 to 191} or {resid 198 to 220}}" frame 0]
    set refprot1 [atomselect 0 "protein and backbone and {{resid 100 to 125}}" frame 0]
    set refprot2 [atomselect 0 "protein and backbone and {{resid 134 to 160}}" frame 0]
    set refprot3 [atomselect 0 "protein and backbone and {{resid 168 to 191}}" frame 0]
    set refprot4 [atomselect 0 "protein and backbone and {{resid 198 to 220}}" frame 0]
    if {$withGBI == 1} {set refgbi [atomselect 0 "segname GBI1 and noh" frame 0] }
    if {$withGBI == 2} {set refgbi [atomselect 0 "segname GBI2 and noh" frame 0] }
    
    # set selections being compared.
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
        [atomselect $molID all] frame $frame
    
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
        lappend listprot $rmsdprot
        lappend listprot1 $rmsdprot1
        lappend listprot2 $rmsdprot2
        lappend listprot3 $rmsdprot3
        lappend listprot4 $rmsdprot4
        if {$withGBI == 1 || $withGBI == 2} {
          set rmsdgbi [measure rmsd $compgbi $refgbi]
          lappend listgbi $rmsdgbi
        }
    }

    # get the average values, excluding reference (RMSD 0)
    if {$withGBI == 1 || $withGBI == 2} {
        set avggbi [average $listgbi]
    } else {
        set avggbi NaN
    }


    set avgprot [average $listprot]
    set avgprot1 [average $listprot1]
    set avgprot2 [average $listprot2]
    set avgprot3 [average $listprot3]
    set avgprot4 [average $listprot4]

    return [list $avgprot $avgprot1 $avgprot2 $avgprot3 $avgprot4 $avggbi]
}



# =============================================================== #

lappend auto_path /data12/cmf/limvt/tempotools/libs
lappend auto_path /home/limvt/Documents/tempotools/libs
package require tempoUserVMD

# load in psf and reference coords
mol addfile $pdb waitfor all
set molID 1  ;# the reference pdb will be molID 0

foreach way [list "F" "R"] {

    # define output file for RMSDs
    set outDataFile [open ${homedir}/02_analysis/2_rmsd/rmsd_avgByWin_byChain-${way}.dat w]
    if {$withGBI == 1 || $withGBI == 2} {
      puts $outDataFile "# Win | helix backbone rmsd | TM helix 1 | TM helix 2 | TM helix 3 | TM helix 4 | 2GBI rmsd (Angstroms)"
    } else {
      puts $outDataFile "# Win | helix backbone rmsd | TM helix 1 | TM helix 2 | TM helix 3 | TM helix 4 (Angstroms)"
    }


    # get list of trajectories, add zero for nums<10
    set dcdlist {}
#    for {set i 1} {$i <= 1} {incr i} 
    for {set i 1} {$i <= $num_wins} {incr i} {
        if {$i < 10} {
          set win 0$i
        } else {
          set win $i
        }
        lappend dcdlist ${homedir}/FEP_${way}/lambda_${win}/alchemy${win}.dcd
    }

    
    # read in each trajectory, calculate RMSD average over all windows, write to file 
    foreach dcd $dcdlist {

        # the main functions
        mol new $psf
        dopbc -file ${dcd} -frames 50:1:249 -mol $molID -ref protein
        set rmsdlist [rmsdTraj $molID $withGBI]

        # write info
        if {$withGBI == 1 || $withGBI == 2} {
          puts $outDataFile "$molID \t [lindex $rmsdlist 0] \t [lindex $rmsdlist 1] \t [lindex $rmsdlist 2] \t [lindex $rmsdlist 3] \t [lindex $rmsdlist 4] \t [lindex $rmsdlist 5]"
        } else {
          puts $outDataFile "$molID \t [lindex $rmsdlist 0] \t [lindex $rmsdlist 1] \t [lindex $rmsdlist 2] \t [lindex $rmsdlist 3] \t [lindex $rmsdlist 4]"
        }

        mol delete $molID
        incr molID
    }

    close $outDataFile



    if {$writeSumTraj} {

        # check that summary trajectory doesn't already exist    
        set outtraj ${homedir}/02_analysis/2_rmsd/lastFrames-${way}.dcd
        if {[file exists $outtraj]} {
            puts "\n\nError, output file already exists: $outtraj\n\n"
            exit
        }

        # write the trajectory
        animate write dcd $outtraj waitfor all
    
        # Record details of summary trajectory in README 
        set fo [open "${homedir}/02_analysis/2_rmsd/README" "a"] 
        puts $fo "\n\nDetails for summary trajectory:"
        puts $fo " - Name:\t$outtraj"
        puts $fo " - Date:\t[clock format [clock seconds]]"
        puts $fo " - Skip:\t$skip"
        puts $fo " - System:\t$pose/$mut/FEP_$way"
        puts $fo " - PSF:\t\t$psf"
        close $fo
    }

}

exit
 
