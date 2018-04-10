

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
    puts "Aligning Hv1 backbone..."
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
    

} ;# end of calc_rmsf


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
    set sel_z0 [split $pre_z0 {,}]
    set sel_z1 [split $pre_z1 {,}]

    # wrap and IGNORE given pdb -- sometimes has error (a=0.000000 b=0.000000 c=0.000000)
    set n [molinfo top get numframes]
    pbc wrap -compound fragment -center com -centersel "protein" -first 1 -last $n ;# zero-based index
    puts "Aligning Hv1 backbone..."
    align_backbone
    puts "Counting waters..."

    # define output file
    set outDataFile [open $outfile w]
    puts $outDataFile "# Frame | number of waters bt Z crds of two sels"
    
    # get reference Z coordinates
    set z0 [[atomselect top $sel_z0 frame 0] get {z}]
    set z1 [[atomselect top $sel_z1 frame 0] get {z}]
    puts $outDataFile "# Selection 1 (z=$z0): $sel_z0\n# Selection 2 (z=$z1): $sel_z1"
    
    # loop over frames
    set frames [molinfo top get numframes]
    set wats [atomselect top "noh and waters and (z<$z1 and z>$z0)"]
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
foreach dcd $dcdlist {
    mol addfile $dcd first 0 last -1 step $inskip waitfor all
}
# wrapping not done here bc rmsd and rmsf use pdb reference, which may not have pbc params

set __after [info procs] ; # get list of avail functions
puts "\n\nAll trajectories loaded; ready for analysis. Available functions:\n[diff $__before $__after]\n\n" 
