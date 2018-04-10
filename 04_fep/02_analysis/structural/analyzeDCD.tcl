

# Purpose: Analyze list of trajectories using RMSD (TODO), RMSF, ...

# Usage: 
#  1. vmd -dispdev none -e file.tcl -args inpsf inskip inpdb indcd [indcd2 indcd3 ...]
#  2. [call some analysis function in this script in VMD terminal/console] 
#
# Full path names for sourcing:
#    gpl: /beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/2_structure/analyzeDCD.tcl
#    cas: /home/limvt/connect/greenplanet/goto-beegfs/hv1/04_fep/analysis/2_structure/analyzeDCD.tcl

# ========================== Variables ========================= #


set inpsf [lindex $argv 0]
set inskip [lindex $argv 1]
set inpdb [lindex $argv 2]

for {set i 3} {$i < $argc} {incr i} {
    lappend dcdlist [lindex $argv $i]
}

# =============================================================== #

package require pbctools
set __before [info procs] ; # get list of avail functions

proc average L {
    # ============================================================
    # Get average of some list of numbers. Used in other functions.
    # ============================================================
    expr ([join $L +])/[llength $L].
}

proc align_backbone {} {
    # ============================================================
    # Align Hv1 protein by backbone of transmembrane regions.
    #   Used in other functions.
    # ============================================================
    set all [atomselect 0 "all"]
    set refprot [atomselect 0 "protein and backbone and {{resid 100 to 125} or {resid 134 to 160} or {resid 168 to 191} or {resid 198 to 220}}" frame 0]
    set compprot [atomselect 0 "protein and backbone and {{resid 100 to 125} or {resid 134 to 160} or {resid 168 to 191} or {resid 198 to 220}}"]
    
    set num_steps [molinfo top get numframes]
    for {set frame 0} {$frame < $num_steps} {incr frame} {
        $compprot frame $frame
        $all frame $frame
        $all move [measure fit $compprot $refprot]
    }
}

proc diff {before after} {
    # ============================================================
    # Extract differences from two lists. Used in this script to
    # list available analysis functions after loading trajectories.
    #
    # References
    #  - https://tinyurl.com/yccerb3k
    # ============================================================
    set result [list]
    foreach name $before {
        set procs($name) 1
    }
    foreach name $after {
        if { ![info exists procs($name)] } {
            lappend result $name
        }
    }
    return [lsort $result]
} ;# end of diff

proc calc_rmsd_hv1 {outfile {level segment} {gbi 0} } {
    # ============================================================
    # Measure RMSD for Hv1. System is aligned by transmembrane
    # backbone, then RMSD can be calculated over backbone, TODO
    # each segment (S1-S4), or average of each residue. TODO
    #
    # Arguments
    #  - outfile : string
    #       Basename of the output files for .dat, .psf, .pdb
    #  - level : string
    #       Level of detail to get RMSD values. Options: backbone, segment, residue
    #  - gbi : integer
    #       2GBI tautomer. 0=absent, 1=taut1, 2=taut2.
    # Returns
    #  - (nothing)
    # Example usage
    #  - calc_rmsd_hv1 rmsd_segmt_f-01 segment 1
    # ============================================================
    global inpsf
    global inpdb
    global dcdlist

    puts "Aligning system by Hv1 transmembrane backbone..."
    align_backbone

    # set groups for RMSD calculation
    puts "Defining groups for RMSD calculation..."
    lappend rgroup "protein and backbone and {{resid 100 to 125} or {resid 134 to 160} or {resid 168 to 191} or {resid 198 to 220}}"
    # determine groups for level of detail for protein
    if {$level == "residue"} {
        for {set resid 88} {$resid < 231} {incr resid} {
            lappend rgroup "protein and backbone and resid $resid"
        }
    } elseif {$level == "segment"} {
        lappend rgroup "protein and backbone and {{resid 100 to 125}}"
        lappend rgroup "protein and backbone and {{resid 134 to 160}}"
        lappend rgroup "protein and backbone and {{resid 168 to 191}}"
        lappend rgroup "protein and backbone and {{resid 198 to 220}}"
    }
    # determine group for ligand if present
    if {$gbi == 1} {
        lappend rgroup "segname GBI1 and noh"
    } elseif {$gbi == 2} {
        lappend rgroup "segname GBI2 and noh"
    }

    # open file for writing output
    set outDataFile [open $outfile.dat w]
    puts $outDataFile "# Data from files:\n#  $inpsf\n#  $dcdlist\n"
    puts $outDataFile "# RMSD (Angstroms)"
    set header "# Win | TM bb"
    if {$level == "residue"} {
        for {set resid 88} {$resid < 231} {incr resid} {
          lappend allres $resid    
        }
        set reslabel [join $allres "\t\t"]
        set header "$header | $reslabel"
    } elseif {$level == "segment"} {
        set header "$header | S1\t\tS2\t\tS3\t\tS4"
    }

    if {$gbi == 1 || $gbi == 2} {
        set header "$header\t\t2GBI"
    }
    puts $outDataFile $header

    # rmsd calculation
    puts "Calculating RMSD..."
    set num_steps [molinfo top get numframes]
    for {set frame 0} {$frame < $num_steps} {incr frame} {
        if {[expr $frame % 50 == 0]} {puts $frame}
        set curr_line "$frame\t"
        foreach x $rgroup {
            set sel0 [atomselect 0 $x frame 0]
            set sel1 [atomselect 0 $x frame $frame]
            set rmsdsel [format "%.4f" [measure rmsd $sel1 $sel0]]
            set curr_line "$curr_line\t$rmsdsel"
        }
        puts $outDataFile $curr_line
    }
    close $outDataFile

} ;# end of calc_rmsd_hv1

proc calc_rmsf_hv1 {outfile} {
    # ============================================================
    # Measure RMSF by residue for Hv1. VMD cannot measure RMSF for
    # group of atoms, so each residue uses dummy atom that is the
    # center of mass of the amino acid (with backbone, no hydrogens)
    #
    # Arguments
    #  - outfile : string
    #       Basename of the output files for .dat, .psf, .pdb
    # Returns
    #  - (nothing)
    # Example usage
    #  - calc_rmsf_hv1 rmsf_hv1
    # References
    #  - rmsf uses average position of frames as reference, https://tinyurl.com/yawa8xjo
    #  - use of dummy atom, https://tinyurl.com/y9l2wkvk
    # ============================================================
    global inpsf
    global inpdb
    global dcdlist

    # open file for writing output
    set outDataFile [open $outfile.dat w]
    puts $outDataFile "# Data from files:\n#  $inpsf\n#  $dcdlist\n"
    puts $outDataFile "# Res | RMSF (Angstroms)"
    
    # rmsf calculation
    set dumid 58230 ; # borrow some lone ion to use as dummy atom
    puts "Aligning system by Hv1 transmembrane backbone..."
    align_backbone
    puts "Calculating RMSF..."
    for {set resid 88} {$resid < 231} {incr resid} {
        set whole [atomselect top "protein and resid $resid"]
        set group [atomselect top "protein and resid $resid and noh"]
        set dummy [atomselect top "index $dumid"]
        
        set num_steps [molinfo top get numframes]
        for {set frame 0} {$frame < $num_steps} {incr frame} {
            $whole frame $frame
            $group frame $frame
            $dummy frame $frame
    
            # measure crds of center of this noh-residue and set crds on dummy atom
            set xyz [measure center $group]
            $dummy set x [lindex $xyz 0]
            $dummy set y [lindex $xyz 1]
            $dummy set z [lindex $xyz 2]
        }
    
        # rmsf calculation over all frames
        set rmsf [measure rmsf $dummy] 
        $whole set occupancy $rmsf
        puts $outDataFile "$resid\t$rmsf"
        
    }

    # write out rmsf info in occupancy column of a PDB file
    animate write pdb $outfile.pdb beg 0 end 0 sel [atomselect top protein]
    animate write psf $outfile.psf beg 0 end 0 sel [atomselect top protein]
    close $outDataFile
    

} ;# end of calc_rmsf_hv1


proc count_wat_z {outfile pre_z0 pre_z1} {
    # ============================================================
    # Count number of waters from -z0 to +z1 coordinate.
    # Specify selection with no spaces, but use commas for multiple words.
    # See example usage.
    #
    # Arguments
    #  - outfile : string
    #       Name of the output file.
    #   - pre_z0 : string
    #       Selection for z0, should be in the lower z plane
    #   - pre_z1 : string
    #       Selection for z1, should be in the upper z plane
    # Returns
    #  - (nothing)
    # Example usage
    #  - count_wat_z waters-in-zrange.dat protein,and,resid,223,and,name,CZ protein,and,resid,208,and,name,CZ
    # References
    #  - http://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/23723.html
    # ============================================================
    set watlist [list]
    set sel_z0 [split $pre_z0 {,}] ;# string selection text
    set sel_z1 [split $pre_z1 {,}]

    # wrap and IGNORE given pdb -- sometimes has error (a=0.000000 b=0.000000 c=0.000000)
    set n [molinfo top get numframes]
    pbc wrap -compound fragment -center com -centersel "protein" -first 1 -last $n ;# zero-based index
    puts "Aligning system by Hv1 transmembrane backbone..."
    align_backbone
    puts "Counting waters..."

    # define output file
    set outDataFile [open $outfile w]
    puts $outDataFile "# Frame | number of waters bt Z crds of two sels"
    
    # get reference Z coordinates
    set z0 [[atomselect top $sel_z0 frame 0] get {z}]
    set z1 [[atomselect top $sel_z1 frame 0] get {z}]
    puts $outDataFile "# Selection 1 (z=$z0): $sel_z0\n# Selection 2 (z=$z1): $sel_z1"

    # get box edges
    set x0 [[atomselect top "protein and resid 173 and name O" frame 0] get {x}]
    set x1 [[atomselect top "protein and resid 134 and name CA" frame 0] get {x}]
    set y0 [[atomselect top "protein and resid 124 and name O" frame 0] get {y}]
    set y1 [[atomselect top "protein and resid 148 and name C" frame 0] get {y}]
    puts $outDataFile "# Edge 1 (x=$x0): $sel_x0\n# Edge 2 (x=$x1): $sel_x1\n# Edge 3 (y=$y0): $sel_y0\n# Edge 4 (y=$y1): $sel_y1"
    
    # loop over frames
    set frames [molinfo top get numframes]
    set wats [atomselect top "noh and waters and (z<$z1 and z>$z0) and (x<$x1 and x>$x0) and (y<$y1 and y>$y0)"]
    for {set i 0} {$i < $frames} {incr i} {
        $wats frame $i
        $wats update
        set num [$wats num]
        lappend watlist $num
        puts $outDataFile "$i $num"
    }
    
    $wats delete
    puts $outDataFile "# --- Average over traj: [average $watlist]"
    close $outDataFile

} ;# end of count_wat_z



# =============================================================== #


# read in data
mol new $inpsf
mol addfile $inpdb        ;# mol 0 == mol top
foreach dcd $dcdlist {    ;# maybe alter the first step to read in if FEP bc 50 frames equil
    mol addfile $dcd first 0 last -1 step $inskip waitfor all
}
# wrapping not done here bc rmsd and rmsf use pdb reference, which may not have pbc params

set __after [info procs] ; # get list of avail functions
puts "\n\nAll trajectories loaded; ready for analysis. Available functions:\n[diff $__before $__after]\n\n" 
