### tcl script to open the GBI configs found in the reverse clustering script.
### Run with VMD in the pdbqt/wt_ans/ directory.

#open file with config pdb names and parse.
set fp [open "configurations.dat" r]
set file_data [read $fp]
set data [split $file_data "\n"]

foreach line $data {
    #get directory and pdb name for each.
    set name [lindex $line 1]
    set dir [lindex [split $name "_"] 0]
    
    #open first pdb as new molec., add the rest to that molec.
    if  {[molinfo top] == -1 } {
	mol new ./${dir}/bs_1/${name}
    } else {
	mol addfile ./${dir}/bs_1/${name}
    }
    mol new ./${dir}/bs_1/${name}   
    unset name
    unset dir
}

close $fp
#now you can visualize the configurations in VMD.
