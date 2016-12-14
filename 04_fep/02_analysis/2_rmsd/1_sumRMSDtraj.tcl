# Purpose of the script: to take and wrap every nth frame for summary trajectory.
# Usage: change files, then "vmdt -e file.tcl"
# vim to replace all foo with bar: :%s/foo/bar/g
# Tips: 
#  - To get endFrames, first need to know number of frames [catdcd -num file.dcd]
#    Then specify frames like [-frames 499:1:500]


# ========================== Variables ========================= #


set pose "17041_19"
set mut "c1_F150A_F182A"
set withGBI 1     ;# 1 for True, 0 for false
set extract last  ;# if you change this, also change "out"
set way "R"
set num_wins 40  
set skip 1


# Edit input psf and output dcd file names.
set dir /data12/cmf/limvt/hv1/04_fep/${pose}/${mut}
set psf ${dir}/00_main/17041_19-F2A.psf                    ;#
set out ${dir}/02_analysis/2_rmsd/finalFrames-${way}.dcd   ;#
set outDataFile [open rmsd_endFrames-${way}.dat w]


# If not 250 frames (5 ns), edit the dopbc line below.



# =============================================================== #
# ============= Grab frame from all windows ===================== #
# =============================================================== #

lappend auto_path /data12/cmf/limvt/tempotools/libs
package require tempoUserVMD

mol new $psf

if {[file exists $out]} {
    puts "\n\nError, output file already exists: $out\n\n"
    exit
}

set dcdlist {}
for {set i 1} {$i <= $num_wins} {incr i} {
    # process numbers < 10 to add a zero
    if {$i < 10} {
        set win 0$i
    } else {
        set win $i
    }

    # get list of trajectories
    lappend dcdlist ${dir}/FEP_${way}/lambda_${win}/alchemy${win}.dcd
}

# add the first frame of the first dcd as reference
dopbc -file [lindex $dcdlist 0] -frames 0:1:0 -ref protein

foreach dcd $dcdlist {
    #dopbc -file ${dcd} -frames 0:$skip:-1 -ref protein ;# get every skipth frame

    if { [string equal -nocase "last" $extract] } {
      dopbc -file ${dcd} -frames 249:1:250 -ref protein ;# get last? frame
    } else {
      dopbc -file ${dcd} -frames 0:1:1 -ref protein     ;# get first frame of each
    }

}

animate write dcd $out waitfor all

set fo [open "${dir}/02_analysis/2_rmsd/README" "a"] 
puts $fo "\n\nDetails for summary trajectory:"
puts $fo " - Name:\t$out"
puts $fo " - Date:\t[clock format [clock seconds]]"
puts $fo " - Skip:\t$skip"
puts $fo " - System:\t$pose/$mut/FEP_$way"
puts $fo " - PSF file:\t$psf"
close $fo

# =============================================================== #
# ============= Evaluate RMSD wrt to frame 1 ==================== #
# =============================================================== #

# file to output data for plotting
if {$withGBI} {
  puts $outDataFile "#Frame | helix backbone rmsd | Res112,150,181,211 rmsd | 2GBI rmsd (Angstroms)"
} else {
  puts $outDataFile "#Frame | helix backbone rmsd | Res112,150,181,211 rmsd (Angstroms)"
}

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
    set rmsdres [measure rmsd $compres $refres]

    if {$withGBI} {
      set rmsdgbi [measure rmsd $compgbi $refgbi]
      puts $outDataFile "$frame \t $rmsdprot \t $rmsdres \t $rmsdgbi"
    } else {
      puts $outDataFile "$frame \t $rmsdprot \t $rmsdres"
    }
}
exit

