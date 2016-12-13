
# View 2GBI in Hv1.
# Usage: "vmd -e file.tcl"


color Display Background white
display projection Orthographic

### load directories
set hdir /work/cluster/limvt/hv1/05_dihedScan/2_tautomer/angles-md
cd $hdir
set dirlist [glob */ ]


### load all coordinates into one mol, like a movie
if {0} {
  mol new ${hdir}/../gbi-tautomer2.psf
  mol delrep 0 0
  mol addrep 0
  mol modselect 0 0 resname GBI2
  mol modstyle 0 0 Licorice 0.200000 12.000000 30.000000
  mol modcolor 0 0 Name
  
  foreach i $dirlist {
      set dir [string trim $i "/"]
      mol addfile ${hdir}/${dir}/gbi2-${dir}.coor
  }

}


### load all coordinates separately to show all frames at once
if {1} {
  set count 0
  foreach i $dirlist {
      set dir [string trim $i "/"]
      mol new ${hdir}/../gbi-tautomer2.psf
      mol addfile ${hdir}/${dir}/gbi2-${dir}.coor
  
      mol delrep 0 $count
      mol addrep $count
      mol modselect 0 $count resname GBI2
      mol modstyle 0 $count Licorice 0.200000 12.000000 30.000000
      mol modcolor 0 $count Name
  
      #mol off $count           ;# don't display anything to start
      incr count 
      cd ../
  }
}

