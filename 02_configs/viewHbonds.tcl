
# View last frame of all taut2 poses
# Arguments: text file with poses, common dcd name, last frame number of common dcd
# vmd -e viewHbonds.tcl -args pose_dirs2.txt npt02.dcd 4999
# vmd -e viewHbonds.tcl

set taut1 1

if {$taut1} {
    cd /pub/limvt/hv1/02_configs/1_tautomer
    set lig GBI1
} else {
    set lig GBI2
}

#color Display Background gray
display projection Orthographic
display depthcue off

# load directories
set fn [lindex $argv 0]

# if no arguments, assume view one molecule.
if { $argc > 0 } {
    set fp [open [lindex $argv 0] r]
    set dirs [read $fp]
    close $fp
} else {
    set dirs "notApplicable"
}


set count 0
foreach pose $dirs {

  if {$argc > 0} {
    cd $pose
    set psf [glob *psf]
    mol new $psf type {psf} first 0 last -1 step 1 waitfor all
    mol addfile [lindex $argv 1] type {dcd} first [lindex $argv 2] last -1 step 1 waitfor -1 $count
  }

    mol delrep 0 $count

    # the ligand
    mol color Name
    mol representation Licorice 0.200000 10.000000 10.000000
    mol selection resname $lig
    mol material Opaque
    mol addrep $count

    # the neighbors
    mol color ResType
    mol representation Licorice 0.100000 10.000000 10.000000
    mol selection (protein or water) and same residue as noh and within 5 of resname $lig
    mol material Opaque
    mol addrep $count

    # the hbond network
    mol color Name
    mol representation HBonds 3.500000 40.000000 6.000000
    mol selection (nitrogen or oxygen) and (resname $lig or (protein or water) and same residue as noh and within 5 of resname $lig)
    mol material Opaque
    mol addrep $count
    
    # phenylalanines
    mol color ColorID 12
    mol representation Licorice 0.200000 10.000000 10.000000
    mol selection protein and resid 149
    mol material Opaque
    mol addrep $count
    # phenylalanines
    mol color ColorID 3
    mol representation Licorice 0.200000 10.000000 10.000000
    mol selection protein and resid 150
    mol material Opaque
    mol addrep $count
    # phenylalanines
    mol color ColorID 6
    mol representation Licorice 0.200000 10.000000 10.000000
    mol selection protein and resid 182
    mol material Opaque
    mol addrep $count

    # hide displays for ease of viewing
    mol off $count           ;# don't display anything to start
    mol showrep $count 1 0   ;# don't show neighbors
    mol showrep $count 2 0   ;# don't show hbonds

    incr count 
    if {$fn != 0 && $taut1} {
        cd ../
    } elseif {$fn != 0} {
        cd ../../
    }
}

display resetview
