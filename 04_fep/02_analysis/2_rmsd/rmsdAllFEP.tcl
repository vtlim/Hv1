

# ====================================================================================
#
# Usage: vmdt -e rmsdAllFEP.tcl -args /work/cluster/limvt/hv1/04_fep/1_tautomer/17041_19/c1_F150A_F182A/
#
# Loop through all frames of all windows for FEP_F and FEP_R and calculate the
#   RMSD of (protein backbone, no loops) and (all specified protein resids).
#   A new *.dat file is created for every window with all specified selections.
#
# Processing later
#   for i in 0{1..9} {10..40}; do echo $i; xg gbi_FEP_F_lambda_$i.dat; done
#
# ====================================================================================


### set list of protein residue numbers (e.g. 211 for R211) to be analyzed. 
### don't use noh because then Hs won't get moved and RMSDs will be inflated
###   when I tested on 19 residues, it took 2 minutes per window. = 160 minutes = almost 2 hours

set prefix selections
#set prefix gbi

set numbers "112 149 150 181 182 185 211 214 113 116 178 207"
#set numbers ""
#set selList ""

#set withLig 0
set withLig 1



# ====================================================================================


### build VMD keywords from residue IDs
foreach id $numbers {lappend selList "protein and resid $id"}
puts $selList

### find system files
set home [lindex $argv 0]
cd $home
set mypsf [glob -dir $home/00_main *psf]
set mypdb [glob -dir $home/00_main *pdb]
set count 1  ;# for mol add,etc. mol 0 is reference.

#### before opening each traj, open pdb reference (!= to frame1 of FEP_F/lambda01)
mol new $mypdb
set refprot [atomselect 0 "protein and backbone and {{resid 100 to 125} or {resid 134 to 160} or {resid 168 to 191} or {resid 198 to 220}}"]
if {$withLig} {set reflig [atomselect 0 "resname GBI2"] }

### go through all fwd and rev trajectories
#foreach group {FEP_R} {
foreach group {FEP_F FEP_R} {

    foreach window [glob $home/$group/lambda_*] {
        set lambda [lindex [split $window /] end]

        ### open file for writing out data
        set outFile [open $home/02_analysis/2_rmsd/${prefix}_${group}_${lambda}.dat w]
        set headers "# prot lig "
        foreach thisSel $numbers {lappend headers "   $thisSel"}
        puts $outFile [join $headers]

        ### open trajectory for this window. not using dopbc since 
        #     only interested in protein RMSD after alignment with pdb. 
        mol new $mypsf
        set dcd [glob -dir $window *dcd]
        mol addfile $dcd type {dcd} first 0 last -1 step 1 waitfor all $count
        set compprot [atomselect $count "protein and backbone and {{resid 100 to 125} or {resid 134 to 160} or {resid 168 to 191} or {resid 198 to 220}}"]
        if {$withLig} {set complig [atomselect $count "resname GBI2"]}
            
        set nframes [molinfo $count get numframes]
        for {set i 0} {$i < $nframes} {incr i} {

            ### calculate tranformation matrix to align protein backbone (no loops) based on $refprot
            $compprot frame $i
            set trans_mat [measure fit $compprot $refprot]

            ### protein: align, compute RMSD
            $compprot move $trans_mat
            set rmsdprot [format "%.2f" [measure rmsd $compprot $refprot]]
            set writeLine "$i $rmsdprot"

            ### ligand: align, compute RMSD
            $complig frame $i
            $complig move $trans_mat
            set rmsdlig [format "%.2f" [measure rmsd $complig $reflig]]
            lappend writeLine " $rmsdlig"


            ### loop over residue selections to compare.
            foreach thisSel $selList {

                # set VMD selections for reference and comparison
                set refsel [atomselect 0 $thisSel]
                set compsel [atomselect $count $thisSel frame $i]

                # selection: align, compute RMSD
                $compsel move $trans_mat
                set rmsdsel [format "%.2f" [measure rmsd $compsel $refsel]]
                lappend writeLine " $rmsdsel"

            } ;# done with this selection
            puts $outFile [join $writeLine]
            unset writeLine


        } ;# done with all frames


        mol delete $count
        incr count

    } ;# done with this window
    close $outFile

} ;# done with this group

exit
