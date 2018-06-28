
# ============================================================================
# After loading psf and dcd, use this script for various visualizations of Hv1.
# Only modifies visualization on TOP molecule.
# Steps:
#   1. "source file.tcl"
#   2. "view_clear"
#   3. "view_[four/protein]"
# By: Victoria Lim
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
    mol modcolor $repcount $moltop ColorID $repcount
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

