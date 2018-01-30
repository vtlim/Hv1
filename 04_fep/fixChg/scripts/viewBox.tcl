# VMD for LINUXAMD64, version 1.9.3 (November 30, 2016)
# Log file 'viewBox.tcl', created by user limvt
mol color Name
mol representation Lines 1.000000
mol selection all
mol material Opaque
mol addrep 0
mol modselect 1 0 protein
mol modstyle 1 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor 1 0 ColorID 0
mol modcolor 1 0 ColorID 11
mol color ColorID 11
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol selection protein
mol material Opaque
mol addrep 0
mol modselect 2 0 protein and resid 124
mol modstyle 2 0 VDW 1.000000 12.000000
mol modcolor 2 0 ColorID 7

source /home/limvt/connect/greenplanet/goto-beegfs/hv1/04_fep/sandbox_chgScript1/drawBox.tcl
draw_box origin 60 60 60
draw_box origin 27 28 35

# VMD for LINUXAMD64, version 1.9.3 (November 30, 2016)
# end of log file.

