# VMD for LINUXAMD64, version 1.9.2 (December 29, 2014)
# Purpose: Check out mutation site of FEP coordinate file.
# Usage: vmd file.fep -e viewMutation.tcl -args PROTresid GBIresname
# Example: vmd file.fep -e viewMutation.tcl -args 150 GBI1

color Display Background white
display projection Orthographic
display depthcue off

set res [lindex $argv 0]
set gbi [lindex $argv 1]
set count 0
mol delrep 0 0

# add representation for file with GBI only
mol addrep $count
mol modselect 0 $count resname $gbi
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

# add representation for overlapping with 2GBI
mol addrep $count
mol modselect 3 $count "not resname $gbi and within 3 of resname $gbi"
mol modstyle 3 $count Licorice 0.300000 12.000000 12.000000
mol modcolor 3 $count ColorID 1

menu graphics on
# VMD for LINUXAMD64, version 1.9.2 (December 29, 2014)
# end of log file.
