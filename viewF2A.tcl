# VMD for LINUXAMD64, version 1.9.2 (December 29, 2014)
# Log file 'viewF150.tcl', created by user victoria

color Display Background white

set res 150

mol delrep 0 0
mol color Name
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol selection protein
mol material Opaque

mol addrep 0
mol modmaterial 0 0 Transparent
mol modcolor 0 0 ColorID 0

mol addrep 0
mol modselect 1 0 resname GBI1
mol modstyle 1 0 Licorice 0.300000 12.000000 12.000000
mol modcolor 1 0 ColorID 0
mol modcolor 1 0 ColorID 7
mol color ColorID 0
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol selection protein
mol material Opaque

mol addrep 0
mol modselect 2 0 protein and resid $res
mol modcolor 2 0 Beta
mol modstyle 2 0 Licorice 0.300000 12.000000 12.000000
mol color ColorID 4
mol representation Licorice 0.300000 12.000000 12.000000
mol selection protein and resid 150
mol material Opaque

# VMD for LINUXAMD64, version 1.9.2 (December 29, 2014)
# end of log file.
