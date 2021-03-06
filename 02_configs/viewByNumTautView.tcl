
# Purpose: View different representations of 2GBI in Hv1.
# Usage: vmd -e file.tcl -args pose_dirs.txt npt01.dcd
# Notes on usage:
#   - The text file should have a list of pose subdirectories inside `tdir`.
#   - The coordinate file is loaded for each of the poses (e.g., npt01 traj, npt03 coor file).
#   - Arguments are optional. Don't include any if you just want to view one mol.
#   - How to view just one mol: (TODO)
#

### Num:  to view one config, use no arguments. Else see example usage.
### Taut: specify set taut1 0 for false (aka taut2) or 1 for true.
### View: options are [four]  for D112,F150,S181,R211 and protein/lig
#                     [hbonds] for H bonding neighborhood
#                     [lipid]  for seeing the membrane/POPC 149
#                     default for the standard black

set taut 2
set tdir /dfs3/pub/limvt/hv1/02_configs/2_tautomer
#set tdir /pub/limvt/hv1/02_configs/1_tautomer
#set tdir /home/limvt/connect/hpc/goto-pub/hv1/02_configs/1_tautomer


set view pose      ;# just 2GBI, colored by pose
#set view four     ;# residues expected to bind 2gbi
#set view hbonds   ;# hbond network close to 2gbi
#set view lipid    ;# carbonyls and hv1
#set view lipid2   ;# lower carbonyls with plane

set frame1 0
set frame2 0
#set frame2 -1
#set frame2 2754
#set frame1 2504
#set frame2 2504

# load files into a single mol
set commonPSF 0
#set commonPSF 1
#set mypsf /pub/limvt/hv1/07_rotateF182/01_setup/hHv1_t1f_npt10.psf
set mypsf /pub/limvt/hv1/07_rotateF182/01_setup/pose17_F150A.psf

# ==========================================

if {$taut == 1} {
    set lig GBI1
} elseif {$taut == 2} {
    set lig GBI2
}
cd $tdir

display projection Orthographic
display depthcue off
menu graphics on

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

    # skip commented out lines
    if {![string match "#*" $pose]} {
        cd $pose
    } else {
        continue
    }


    # load files into a new mol each
    if {! $commonPSF} {
        set mypsf [glob *psf]
        mol new $mypsf type {psf} first 0 last -1 step 1 waitfor all
    } elseif {$commonPSF && $count == 0} {
        mol new $mypsf type {psf} first 0 last -1 step 1 waitfor all
        set origRep 1
    } else {
        set count 0
        set origRep 0
    }


    mol addfile [lindex $argv 1] first $frame1 last $frame2 step 1 waitfor -1 $count
    #mol addfile [lindex $argv 1] type {dcd} first $frame1 last $frame2 step 1 waitfor -1 $count
    #mol addfile [glob *wGBI.pdb]
    mol rename $count $pose
  }

    if {($view == "hbonds" && !$commonPSF) || ($view == "hbonds" && $origRep == 1)} {

        mol delrep 0 $count
        color Display Background gray

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
        #mol off $count           ;# don't display anything to start
        #mol showrep $count 1 0   ;# don't show neighbors
        #mol showrep $count 2 0   ;# don't show hbonds
        mol showrep $count 3 0   ;# don't show Phe
        mol showrep $count 4 0   ;# don't show
        mol showrep $count 5 0   ;# don't show
        display resetview

    } elseif {$view == "lipid"} {

        # 1 protein secondary structure
        mol color Structure
        mol representation NewCartoon 0.300000 10.000000 4.100000 0
        mol selection protein
        mol material Opaque
        mol addrep $count

        # 2 lipid carbonyl carbons
        mol color Structure
        mol representation VDW 1.000000 12.000000
        mol selection name C21 C31
        mol material Opaque
        mol addrep $count

        # 3 rogue lipid
        mol color Name
        mol representation VDW 1.000000 12.000000
        mol selection lipid and resid 149
        mol material Opaque
        mol addrep $count

        # 4 lipids close to 2gbi
        mol color Structure
        mol representation VDW 1.000000 12.000000
        mol selection lipid and within 5 of resname $lig
        mol material Opaque
        mol addrep $count

        # 5 2gbi itself
        mol color Name
        mol representation Licorice 0.300000 10.000000 10.000000
        mol selection resname $lig
        mol material Opaque
        mol addrep $count


        # 6 add representation for D112
        mol color ColorID 10
        mol representation Licorice 0.300000 10.000000 10.000000
        mol selection protein and resid 112
        mol material Opaque
        mol addrep $count

        # 7 add representation for F150
        mol color ColorID 3
        mol representation Licorice 0.300000 10.000000 10.000000
        mol selection protein and resid 150
        mol material Opaque
        mol addrep $count

        # 8 add representation for S181
        mol color ColorID 7
        mol representation Licorice 0.300000 10.000000 10.000000
        mol selection protein and resid 181
        mol material Opaque
        mol addrep $count

        # 9 add representation for R211
        mol color ColorID 11
        mol representation Licorice 0.300000 10.000000 10.000000
        mol selection protein and resid 211
        mol material Opaque
        mol addrep $count

        # 10 F149
        mol color ColorID 12
        mol representation Licorice 0.200000 10.000000 10.000000
        mol selection protein and resid 149
        mol material Opaque
        mol addrep $count

        # 11 F182
        mol color ColorID 6
        mol representation Licorice 0.200000 10.000000 10.000000
        mol selection protein and resid 182
        mol material Opaque
        mol addrep $count

        mol off $count    ;# don't display anything to start
        mol showrep $count 0 0 ;# for countth mol, 0th rep, turn off view
        mol showrep $count 5 0 ;# hide 2gbi
        mol showrep $count 6 0
        mol showrep $count 7 0
        mol showrep $count 8 0
        mol showrep $count 9 0
        mol showrep $count 10 0
        mol showrep $count 11 0
        display resetview

    } elseif {$view == "lipid2"} {

        # 0 protein secondary structure
        mol color Structure
        mol representation NewCartoon 0.300000 10.000000 4.100000 0
        mol selection protein
        mol material Opaque
        mol addrep $count

        # 1 rogue lipid
        mol color Name
        mol representation VDW 1.000000 12.000000
        mol selection lipid and resid 149
        mol material Opaque
        mol addrep $count

        # 2 lipid carbonyl carbons
        mol color Structure
        mol representation VDW 1.000000 12.000000
        mol selection name C21 C31 and z < 0
        mol material Opaque
        mol addrep $count

        # 3 intrusive carbonyl carbons
        mol color ColorID 3
        mol representation VDW 1.000000 12.000000
        mol selection resid 149 and (name C21 or name C31)
        mol material Opaque
        mol addrep $count

        # 4 add representation for N214
        mol color ColorID 7
        mol representation Licorice 0.300000 10.000000 10.000000
        mol selection protein and resid 214
        mol material Opaque
        mol addrep $count

        # 5 add representation for R211
        mol color ColorID 11
        mol representation Licorice 0.300000 10.000000 10.000000
        mol selection protein and resid 211
        mol material Opaque
        mol addrep $count

        # 6 R226
        mol color ColorID 11
        mol representation Licorice 0.300000 10.000000 10.000000
        mol selection protein and resid 226
        mol material Opaque
        mol addrep $count


        #mol off $count    ;# don't display anything to start
        mol showrep $count 0 0 ;# for countth mol, 0th rep, turn off view
        mol showrep $count 1 0 ;# hide 2gbi
        mol showrep $count 5 0
        mol showrep $count 6 0

        # draw plane for reference (2 triangles)
        draw color green
        draw triangle {-50 -40 -8} {-50 50 -8} {50 -40 -8}
        draw triangle {-50 50 -8} {50 50 -8} {50 -40 -8}
        display resetview

    } elseif {$view == "four"} {


        color Display Background white

        # add representation for file with GBI only
        mol modselect 0 $count resname $lig
        mol modstyle 0 $count Licorice
        mol modcolor 0 $count ColorID 4
        #mol modcolor 0 $count Name

        # add representation for protein
        mol addrep $count
        mol modselect 1 $count protein
        mol modstyle 1 $count NewCartoon
        mol modcolor 1 $count ColorID 0
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

        # add representation for residues
        mol addrep $count
        mol modselect 6 $count protein and resid 182
        mol modstyle 6 $count Licorice 0.300000 12.000000 12.000000
        mol modcolor 6 $count ColorID 6

        # add representation for residues
        mol addrep $count
        mol modselect 7 $count protein and resid 185
        mol modstyle 7 $count Licorice 0.300000 12.000000 12.000000
        mol modcolor 7 $count ColorID 23

        # add representation for overlapping with 2GBI
        mol addrep $count
        mol modselect 8 $count "not resname $lig and within 3 of resname $lig"
        mol modstyle 8 $count Licorice 0.300000 12.000000 12.000000
        mol modcolor 8 $count ColorID 1

        # rogue lipid
        mol addrep $count
        mol modselect 9 $count "lipid and resid 149"
        mol modstyle 9 $count VDW 1.000000 12.000000
        mol modcolor 9 $count Name

        # hide displays for ease of viewing
        #mol off $count           ;# don't display anything to start
        mol showrep $count 1 0   ;# don't show protein
        mol showrep $count 6 0   ;# don't show resid 182
        mol showrep $count 7 0   ;# don't show resid 185
        mol showrep $count 8 0   ;# don't show 2GBI close environment
        mol showrep $count 9 0   ;# don't show lipid 149
        display resetview

    } elseif {$view == "pose"} {

        color Display Background white

        # add representation for file with GBI only
        mol modselect 0 $count resname $lig
        mol modstyle 0 $count Licorice
        mol modcolor 0 $count Molecule

        # add representation for protein
        mol addrep $count
        mol modselect 1 $count protein
        mol modstyle 1 $count NewCartoon
        mol modcolor 1 $count ColorID 0
        mol modmaterial 1 $count Transparent

        # add representation for overlapping with 2GBI
        mol addrep $count
        mol modselect 2 $count "not resname $lig and within 3 of resname $lig"
        mol modstyle 2 $count Licorice 0.300000 12.000000 12.000000
        mol modcolor 2 $count ColorID 1

        # hide displays for ease of viewing
        #mol off $count           ;# don't display anything to start
        mol showrep $count 1 0   ;# don't show protein
        mol showrep $count 2 0   ;# don't show 2GBI close environment
        display resetview

    }


    incr count

    if {$taut == 1} {
        cd ../
    } elseif {$taut == 2} {
        cd ../../
    }
}


# =====================
# let's say you add a set of npt03.coor's and want to add additional files to each
# =====================
# set count 0
# set fp [open pose_dirs.txt r]
# set dirs [read $fp]
# close $fp
#
# foreach pose $dirs {
#   cd $pose
#   mol addfile hHv1_open_wGBI.pdb $count
#   cd ../
#   incr count
# }
#
