### tcl script to perform reverse clustering on all of the 2GBI docking poses 
### in monomer open-state Hv1. Look through 400 poses and cluster them based 
### on RMSD. Then choose 1 out of each cluster to use in subsequent MD sim 
### of Hv1+GBI. Run with VMD in pdbqt/wt_ans/ directory.

#open all 400 docking poses in VMD. The Hv1 configs were aligned prior to
#the docking simulation, so don't need to align 2GBI poses.
set dir [glob -type d hv1-fit-*]

#loop through the directory containing the docked poses
foreach pose $dir {
    
    #open all 20 poses per protein config
    for {set i 1} {$i < 21} {incr i} {
	set num [format "%02d" $i]
	
	#open first pdb as new molec, add the rest to that molec.
	if {[molinfo top] == -1 } {
	    mol new ./${pose}/bs_1/${pose}_bs1_docked_ligand_${num}.pdb
	} else {
	    mol addfile ./${pose}/bs_1/${pose}_bs1_docked_ligand_${num}.pdb
	}
	
	#make a list of the filenames to match to frame#
	lappend filelist "${pose}_bs1_docked_ligand_${num}.pdb"
	unset num
    }
    puts "Loaded pose ${pose} ..."
}

#use VMD measure cluster command to cluster the 400 configs based on RMSD
set sel [atomselect top "noh"]
set cluster [measure cluster $sel first 1 last -1 distfunc rmsd num 50 cutoff 3.0]
 
#find 1 config from each cluster to do quick MD sim with.
#write the config # to file.
set fp [open "configurations.dat" "w"]
for {set j 0} {$j < 51} {incr j} {
    set conNum [lindex [lindex $cluster $j] 0]
    set conID [lindex $filelist [expr $conNum]] 
    puts $fp "Cluster${j}  ${conID}"
}
close $fp


exit
