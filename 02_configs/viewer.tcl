
# ============================================================================
# After loading psf and dcd, use this script for various visualizations of Hv1.
# Only modifies visualization on TOP molecule.
# Steps:
#   1. "source file.tcl"
#   2. "view_clear"
#   3. "view_[four/protein]"
# By: Victoria Lim
#
# Source:
#   * Greenplanet: /beegfs/DATA/mobley/limvt/hv1/github/02_configs/viewer.tcl
#   * Cassandra:   /home/limvt/connect/greenplanet/goto-beegfs/hv1/github/02_configs/viewer.tcl
# ============================================================================

#set psf [glob *psf]
#mol new $psf type {psf} first 0 last -1 step 1 waitfor all
#mol addfile npt_eq01.dcd type {dcd} first 0 last -1 step 1 waitfor -1 $moltop


# VMD display settings
color Display Background white
display projection Orthographic
display depthcue off
axes location LowerLeft
set repcount 1

proc view_clear {} {
    foreach moltop [molinfo list] {
        # assumes no more than 20 representations drawn
        for {set i 0} {$i < 20} {incr i} {
            mol delrep 0 $moltop
        }
    }
    set ::repcount 0
}

proc view_protein {} {
    global repcount
    set moltop [molinfo top]

    # add representation for file with GBI only
    mol addrep $moltop
    mol modselect $repcount $moltop resname GBI1
    mol modstyle $repcount $moltop Licorice
    mol modcolor $repcount $moltop ColorID 4
    incr ::repcount

    # add representation for protein
    mol addrep $moltop
    mol modselect $repcount $moltop protein
    mol modstyle $repcount $moltop NewCartoon
    mol modcolor $repcount $moltop ColorID 6
    #mol modmaterial 1 $moltop Transparent
    incr ::repcount

    # add representation for S1
    mol addrep $moltop
    mol modselect $repcount $moltop protein and resid 99 to 125
    mol modstyle $repcount $moltop NewCartoon
    mol modcolor $repcount $moltop ColorID 10
    incr ::repcount

    # add representation for S2
    mol addrep $moltop
    mol modselect $repcount $moltop protein and resid 134 to 160
    mol modstyle $repcount $moltop NewCartoon
    mol modcolor $repcount $moltop ColorID 3
    incr ::repcount

    # add representation for S3
    mol addrep $moltop
    mol modselect $repcount $moltop protein and resid 168 to 191
    mol modstyle $repcount $moltop NewCartoon
    mol modcolor $repcount $moltop ColorID 7
    incr ::repcount

    # add representation for S4
    mol addrep $moltop
    mol modselect $repcount $moltop protein and resid 198 to 222
    mol modstyle $repcount $moltop NewCartoon
    mol modcolor $repcount $moltop ColorID 11
    incr ::repcount

}

proc view_four {} {
    global repcount
    set moltop [molinfo top]

    # add representation for file with GBI only
    mol addrep $moltop
    mol modselect $repcount $moltop resname GBI1
    mol modstyle $repcount $moltop Licorice
    mol modcolor $repcount $moltop ColorID 4
    incr ::repcount

    # add represe$repcount ntation for protein
    mol addrep $moltop
    mol modselect $repcount $moltop protein
    mol modstyle $repcount $moltop NewCartoon
    mol modcolor $repcount $moltop ColorID 0
    mol modmaterial $repcount $moltop Transparent
    incr ::repcount

    # add represe$repcount ntation for residues
    mol addrep $moltop
    mol modselect $repcount $moltop protein and resid 112
    mol modstyle $repcount $moltop Licorice 0.300000 12.000000 12.000000
    mol modcolor $repcount $moltop ColorID 10
    incr ::repcount

    # add represe$repcount ntation for residues
    mol addrep $moltop
    mol modselect $repcount $moltop protein and resid 150
    mol modstyle $repcount $moltop Licorice 0.300000 12.000000 12.000000
    mol modcolor $repcount $moltop ColorID 3
    incr ::repcount

    # add represe$repcount ntation for residues
    mol addrep $moltop
    mol modselect $repcount $moltop protein and resid 181
    mol modstyle $repcount $moltop Licorice 0.300000 12.000000 12.000000
    mol modcolor $repcount $moltop ColorID 7
    incr ::repcount

    # add represe$repcount ntation for residues
    mol addrep $moltop
    mol modselect $repcount $moltop protein and resid 211
    mol modstyle $repcount $moltop Licorice 0.300000 12.000000 12.000000
    mol modcolor $repcount $moltop ColorID 11
    incr ::repcount

    # hide displays for ease of viewing
    #mol off $moltop           ;# don't display anything to start
    #mol showrep $moltop 1 0   ;# don't show protein
}

proc view_dens_wat { {infile "watdens.dx"} } {
    # ============================================================
    # View volumetric density of water in given DX file.
    #
    # Arguments
    #  - infile : string
    #      Name of the output file. Default is "watdens.dx".
    # Returns
    #  - (nothing)
    # Example usage
    #  - calc_dens_wat watdens.dx
    #  - calc_dens_wat
    # Notes
    #  - If you want to use this function, the original call to analyzeDCD should not be vmdt.
    # ============================================================
    mol new $infile type {dx} first 0 last -1 step 1 waitfor 1 volsets {0 }
    set moltop [molinfo top]

    mol addrep $moltop
    mol modstyle 0 $moltop Isosurface
    mol modcolor 0 $moltop Volume 0
    mol modstyle 0 $moltop Isosurface 0.032205 0 2 1 1 1
    mol modstyle 0 $moltop Isosurface 0.032205 0 0 1 1 1
    mol scaleminmax $moltop 0 -0.150000 0.050000 ;# note that moltop and repnum are switched

} ;# end of view_dens_wat
