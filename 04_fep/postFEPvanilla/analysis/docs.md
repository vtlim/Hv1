/* vim: set filetype=markdown : */


# Documentation for analysis scripts for FEP simulations in NAMD
Last updated: Apr 10 2018

-----------------------------------------------------------------------------------------------------------------

## Contents

* `2_rmsd`
* `3_ligClose`
* `4_nonCovInts`
* `5_rmsf`
* `6_dists`
* `7_waters`

-----------------------------------------------------------------------------------------------------------------

Documentation template:

## `Analysis directory`

Filename(s) go here.
* *Description*:
* *Parameters*:
* *Returns*:
* *Usage*:
* *Example*:

-----------------------------------------------------------------------------------------------------------------


## View scripts

`/beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/viewHbonds.vmd`
* *Description*:
  - View hbond contacts for some selection.
* *Example*:
  - `vmd -e /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/viewHbonds.vmd -args ../taut1_18629-19_k208restr.psf ../npt01-fwd.dcd`


-----------------------------------------------------------------------------------------------------------------


## `2_rmsd`
`/beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/structural/analyzeDCD.tcl`  
`/DFS-L/DATA/mobley/limvt/analysis/plotXY.py`

* *Description*:
   * Calculate RMSD of Hv1 along trajectory, by Hv1 helix and overall (includes 2GBI if present)
* *Examples*:
   * `vmdt -e /beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/structural/analyzeDCD.tcl -args ../../taut1_18629-19_s211.psf 1 ../../taut1_18629-19_s211.pdb ../../npt01-fwd.dcd`
      - `calc_rmsd_hv1 rmsd_hv1 segment 1`
      - `calc_rmsf_hv1 rmsf_hv1`
   * `python /DFS-L/DATA/mobley/limvt/analysis/plotXY.py -i rmsd_hv1.dat -m 100 -x 'time (ns)' -y 'RMSD (A)' -l 'TM bb;Hv1 S1;Hv1 S2;Hv1 S3;Hv1 S4;2GBI (t1)' -o rmsd_chain-runAvg.png`
   * TODO `python /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/1_rmsd/xyPlot.py -i rmsd_byChain.dat -c 1,6 -m 100 -x 'time (ns)' -y 'RMSD (A)' -l 'Hv1 backbone;2GBI (t1)' -o rmsd_gbi-runAvg.png`


-----------------------------------------------------------------------------------------------------------------


## `3_ligClose`
* *Description*:
   * `lig2protDist.py`
   * `ligCloseWat.tcl`
   * `ligHbonds.tcl`
* *Usage*:
   * `vmdt -e /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/3_ligClose/ligCloseWat.tcl -args psf dcd sel`
* *Examples*:
   * `vmdt -e /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/3_ligClose/ligCloseWat.tcl -args ../../*.psf ../../npt01.dcd resname,GBI2`
   * `vmdt -e /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/3_ligClose/ligHbonds.tcl -args ../../*.psf ../../npt01.dcd resname,GBI1`


-----------------------------------------------------------------------------------------------------------------


## `4_nonCovInts`

* *Description*:
   * Analyze contact interactions between Hv1/2GBI system.
   * `getEdges.v5.tcl`
   * `getNodes.v5.tcl`
   * `nodeContacts.py`
   * `node.slurm`
* *Notes*:
   * See README inside directory for usage.
   * Original documentation: `/beegfs/DATA/mobley/limvt/hv1/04_fep/compare_r208k/nodeAnalysis/README.md`
   * See subgroup data from 2017-09-14 for R208K.
   * See past applications:
      * `/beegfs/DATA/mobley/limvt/hv1/04_fep/compare_r208k/nodeAnalysis/`
      * `/beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/a1_D112E/1_E112/analysis/4_nonCovInts/`
      * `/beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/b1_V178A/1_A178/analysis/4_nonCovInts/` 


-----------------------------------------------------------------------------------------------------------------


## `5_rmsf`
`/beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/structural/analyzeDCD.tcl`  
`/DFS-L/DATA/mobley/limvt/analysis/plotXY.py`

* *Description*:
   * Calculate RMSF of Hv1 residues.
* *Examples*:
   * `vmdt -e /beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/structural/analyzeDCD.tcl -args ../../taut1_18629-19_s211.psf 1 ../../taut1_18629-19_s211.pdb ../../npt01-fwd.dcd`
      - `calc_rmsf_hv1 rmsf_hv1`
   * `python /DFS-L/DATA/mobley/limvt/analysis/plotXY.py -i rmsf_hv1.dat -m 10`


-----------------------------------------------------------------------------------------------------------------


## `6_dists`

`/beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/6_dists/calcDist.tcl`
* *Description*: 
   * Calculate the distance between one selection of interest and a number of other selections of interest.
* *Parameters*: 
   * PSF file
   * DCD file
   * main selection of interest, 
   * other selection to measure distance from. 
* *Returns*: 
   * Text file with info on input parameters, and columns of (1) frame number, (2) distance1, (3) distance2, etc.  
     Output file is named measDist.dat so rename if rerunning the script on new input data.
* *Usage*: 
   * `vmd -dispdev none -e calcDist.tcl -args input.psf input.dcd atom-selection1 atom-selection2 [atom-selection3] [atom-selection4]`   
* *Examples*:
   * `vmdt -e calcDist.tcl -args hHv1_open_wGBI.psf npt10.dcd protein,and,resid,211,and,name,CZ protein,and,resid,150,and,name,CG`
* *Notes*: 
   * If there are spaces in the selection, separate them with COMMAS, else will be parsed as multiple selections.
     (Quotation marks, square brackets, and curly brackets don't help.)


`/beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/6_dists/calcAng.tcl`
* *Description*: 
   * Measure dihedral angle of atom1-atom2-atom3-atom4.
* *Parameters*: 
   * PSF file
   * DCD file
   * main selection of interest, 
   * other selection to measure distance from. 
* *Returns*: 
* *Usage*: 
* *Examples*:
* *Notes*: 
   - See note for calcDist.tcl.


`/beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/6_dists/polarPlot.py`
* *Description*: 
   * Generate scatter plot or histogram over a polar coordinate plot. Ideal for angle measurement data.
* *Usage*: 
* *Examples*:

-----------------------------------------------------------------------------------------------------------------


## `7_waters`


-----------------------------------------------------------------------------------------------------------------
