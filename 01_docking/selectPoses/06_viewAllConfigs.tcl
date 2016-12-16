
# View 2GBI in Hv1.
# Usage: "vmd -e file.tcl"

cd z_togetherCentered

set dirs [glob *pdb]

color Display Background white
display projection Orthographic

# open and set representations for each
set count 0

puts [lindex $dirs 0]
foreach pose $dirs {
    puts $pose
    mol new $pose first 0 last -1 step 1 waitfor all
    mol delrep 0 $count
    
    # add representation for GBI
    mol addrep $count
    mol modselect 0 $count resname GBI
    mol modstyle 0 $count Licorice
    #mol modcolor 0 $count ColorID [expr $count % 33]
    mol modcolor 0 $count ColorID 4

    # add representation for protein
    mol addrep $count
    mol modselect 1 $count protein
    mol modstyle 1 $count NewCartoon
    mol modcolor 1 $count ColorID [expr $count % 33]
    mol modmaterial 1 $count Transparent

    # add representation for residues
    mol addrep $count
    mol modselect 2 $count protein and resid 112
    mol modstyle 2 $count Licorice 0.300000 12.000000 12.000000
    mol modcolor 2 $count ColorID 10

    # add representation for residues
    mol addrep $count
    mol modselect 3 $count protein and resid 150
    mol modstyle 3 $count Licorice 0.300000 12.000000 12.000000
    mol modcolor 3 $count ColorID 3

    # add representation for residues
    mol addrep $count
    mol modselect 4 $count protein and resid 181
    mol modstyle 4 $count Licorice 0.300000 12.000000 12.000000
    mol modcolor 4 $count ColorID 7

    # add representation for residues
    mol addrep $count
    mol modselect 5 $count protein and resid 211
    mol modstyle 5 $count Licorice 0.300000 12.000000 12.000000
    mol modcolor 5 $count ColorID 11


    # hide displays for ease of viewing
    mol off $count           ;# don't display anything to start
    #mol showrep $count 1 0   ;# don't show protein
    incr count 
}






