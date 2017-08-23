
# Purpose: Compute RMSD of Hv1 model between two coordinates, such as before
#     and after FEP mutation. Protein is initially aligned via transmembrane
#     helices (skipping loops).
# Usage: 
#   - vmdt -e file.tcl -args pdb1 pdb2 beg end skip
#   - beg is start residue number,
#   - end is end residue number
#   - skip is single residue number to not include (as in case of mutation)
#
# To use as a standalone function with command line:
#   - open structures in VMD
#   - source path/with/script.tcl
#   - compare2RMSD molID0 molid1 begresid endresid skipID
#   - compare2RMSD 0 1 88 230 178
#
# TODO: 2GBI RMSD
#   - withGBI: 0 noGBI, 1 taut1, 2 taut2
#
# Author: Victoria Lim | Mobley Lab @ UCI



# ========================== Variables ========================= #


set pdb0 [lindex $argv 0]
set pdb1 [lindex $argv 1]
set beg [lindex $argv 2]
set end [lindex $argv 3]
set skipID [lindex $argv 4]


# ============================================================== #


## specify ligand
#if {$withGBI == 1} {
#    set GBI "resname GBI1 and noh"
#} elseif {$withGBI == 2} {
#    set GBI "resname GBI2 and noh"
#} else {
#    puts "Specify a valid version of 2GBI ligand."
#    exit
#}

# open the trajectory to run the RMSD calculation.
if { $pdb0 ne "" } {
    mol new $pdb0
} else {
    puts "No variable specified for pdb0"
}
if { $pdb1 ne "" } {
    mol new $pdb1
} else {
    puts "No variable specified for pdb1"
}

proc compare2RMSD {molid0 molid1 beg end skipID} {
    set align "protein and backbone and not resid $skipID and not {{resid 126 to 133} or {resid 161 to 167} or {resid 192 to 197} or {resid 221 to 230}}"
    
    # file to output data for plotting
    set outDataFile [open residrmsd.dat w]
    puts $outDataFile "# resid | rmsd (Angs) | resname"
    
    # align the system bring pdb1 protein to pdb0 protein
    set all0 [atomselect $molid0 all]
    set all1 [atomselect $molid1 all]
    set ref0 [atomselect $molid0 "$align"]
    set ref1 [atomselect $molid1 "$align"]
    set trans_mat [measure fit $ref0 $ref1]
    $all0 move $trans_mat  ;# tested interactively
    
    # calc rmsd for each residue
    for {set resid $beg} {$resid < $end} {incr resid} {
        if {$resid == $skipID} {
            continue
        }
        set sel0 [atomselect $molid0 "protein and resid $resid and noh"]
        set sel1 [atomselect $molid1 "protein and resid $resid and noh"]
        set irmsd [measure rmsd $sel0 $sel1]

        # get resname for printing
        set rname [lindex [$sel0 get resname] 0]
        puts $outDataFile "$resid\t$irmsd\t#$rname"
    }
    
    close $outDataFile
}

if { $pdb1 ne "" } {
    compare2RMSD 0 1 $beg $end $skipID
    exit
}
