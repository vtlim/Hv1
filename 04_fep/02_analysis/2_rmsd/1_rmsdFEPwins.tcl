
# Purpose: measure backbone RMSD of every last (or otherwise specified) frame 
#          of all fwd and rev FEP windows compared to initial seed coordinates (PDB).

# Usage: vmdt -e file.tcl -args homedir withGBI
#   - homedir: full path ending with pose and mutation, e.g. /path/to/19415-10/a1_F150A
#   - withGBI: 0 False, 1 True
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
set fo [open "${homedir}/02_analysis/2_rmsd/README" "a"] 



# Archived.
#set skip 1



# =============================================================== #
# ============= Grab frame from all windows ===================== #
# =============================================================== #

lappend auto_path /data12/cmf/limvt/tempotools/libs
lappend auto_path /home/limvt/Documents/tempotools/libs

package require tempoUserVMD
set count 0

# loop over forward and reverse trajectories
foreach way [list "F" "R"] {

    # load in psf and reference coords
    mol new $psf
    mol addfile $pdb

    # define output files
    set outtraj ${homedir}/02_analysis/2_rmsd/lastFrames-${way}.dcd   ;#
    set outDataFile [open ${homedir}/02_analysis/2_rmsd/rmsd_lastFrames-${way}.dat w]

    # check that summary trajectory doesn't already exist    
#    if {[file exists $outtraj]} {
#        puts "\n\nError, output file already exists: $outtraj\n\n"
#        exit
#    }
    
    set dcdlist {}
    for {set i 1} {$i <= $num_wins} {incr i} {

        # process numbers < 10 to add a zero
        if {$i < 10} {
            set win 0$i
        } else {
            set win $i
        }
    
        # get list of trajectories
        lappend dcdlist ${homedir}/FEP_${way}/lambda_${win}/alchemy${win}.dcd
    }
    
    
    foreach dcd $dcdlist {
        #dopbc -file ${dcd} -frames 0:$skip:-1 -ref protein ;# get every skipth frame
        mol addfile ${dcd} first 250 last 250 step 1 waitfor -1 $count
        dopbc -file ${dcd} -frames 249:1:250 -ref protein ;# get last frame (assuming 250 total)
    
    }
    
    # Write out summary trajectory - optional, for later viewing purposes
    animate write dcd $outtraj waitfor all

    # Record details of summary trajectory in README 
    puts $fo "\n\nDetails for summary trajectory:"
    puts $fo " - Name:\t$outtraj"
    puts $fo " - Date:\t[clock format [clock seconds]]"
    #puts $fo " - Skip:\t$skip"
    puts $fo " - System:\t$pose/$mut/FEP_$way"
    puts $fo " - PSF:\t\t$psf"
    
    # =============================================================== #
    # ============= Evaluate RMSD wrt to frame 1 ==================== #
    # =============================================================== #
    
    # file to output data for plotting
    if {$withGBI} {
      puts $outDataFile "# Frame | helix backbone rmsd | 2GBI rmsd (Angstroms)"
    } else {
      puts $outDataFile "# Frame | helix backbone rmsd (Angstroms)"
    }
    
    # set frame 0 as the reference
    set refprot [atomselect top "protein and backbone and {{resid 100 to 125} or {resid 134 to 160} or {resid 168 to 191} or {resid 198 to 220}}" frame 0]
    set refgbi [atomselect top "segname GBI1 and noh" frame 0]
    
    # selections being compared.
    set compprot [atomselect top "protein and backbone and {{resid 100 to 125} or {resid 134 to 160} or {resid 168 to 191} or {resid 198 to 220}}"]
    set compgbi [atomselect top "segname GBI1 and noh"]
    
    set all [atomselect top all]
    
    # calc rmsd for each frame.
    set num_steps [molinfo $count get numframes]
    for {set frame 0} {$frame < $num_steps} {incr frame} {
    
        # get frame of interest
        $compprot frame $frame
        $compgbi frame $frame
    
        # compute transformation
        set trans_mat [measure fit $compprot $refprot]
    
        # do the alignment
        $compprot move $trans_mat
        $compgbi move $trans_mat
    
        # compute the RMSD
        set rmsdprot [measure rmsd $compprot $refprot]
    
        if {$withGBI} {
          set rmsdgbi [measure rmsd $compgbi $refgbi]
          puts $outDataFile "$frame \t $rmsdprot \t \t $rmsdgbi"
        } else {
          puts $outDataFile "$frame \t $rmsdprot \t"
        }
    }
    mol delete 0
    incr count
}

close $fo
#exit
 
