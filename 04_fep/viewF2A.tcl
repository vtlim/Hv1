# VMD for LINUXAMD64, version 1.9.2 (December 29, 2014)
# Log file 'viewF150.tcl', created by user victoria

color Display Background white
display projection Orthographic

set res 150
set count 0
mol delrep 0 0

# add representation for file with GBI only
mol addrep $count
mol modselect 0 $count resname GBI1
mol modstyle 0 $count Licorice
mol modcolor 0 $count ColorID 4

# add representation for protein
mol addrep $count
mol modselect 1 $count protein
mol modstyle 1 $count NewCartoon
mol modcolor 1 $count ColorID 0
mol modmaterial 1 $count Transparent

# add representation for residue
mol addrep $count
mol modselect 2 $count protein and resid $res
mol modstyle 2 $count Licorice 0.300000 12.000000 12.000000
mol modcolor 2 $count Beta

# VMD for LINUXAMD64, version 1.9.2 (December 29, 2014)
# end of log file.
