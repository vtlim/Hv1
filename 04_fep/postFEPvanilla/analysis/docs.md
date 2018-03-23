/* vim: set filetype=markdown : */


# Documentation for analysis scripts for FEP simulations in NAMD
Last updated: Mar 23 2018

Documentation template:

### `title.ext`
* *Description*:
* *Parameters*:
* *Returns*:
* *Usage*:
* *Example*:

-----------------------------------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------------------------------
*View contacts to some selection:*
```
#vmd -e /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/viewHbonds.vmd -args ../taut1_18629-19_k208restr.psf ../npt01-fwd.dcd
```


### `1_rmsd`
`/beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/1_rmsd/rmsd_trajectory.tcl`
`/beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/1_rmsd/xyPlot.py`
`/beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/2_rmsd/2_plotTrajDat.py` ... xyPlot.py may be a potential redundancy of this script??

* *Description*:
   * Calculate RMSD of Hv1 along trajectory, by Hv1 helix and overall (includes 2GBI if present)
   * `rmsd_trajectory.tcl`
   * `xyPlot.py`
* *Usage*:
   * `vmdt -e rmsd_trajectory.tcl -args pdbref withGBI inpsf indcd`
* *Examples*:
   * Easier to copy and modify the runcmd document for each case. Change file names as well as GBI type.
   * `cd /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/f1_R211S/1_S211/analysis/1_rmsd`
   * `vmdt -e rmsd_trajectory.tcl -args ../../taut1_18629-19_s211.pdb 1 ../../taut1_18629-19_s211.psf ../../npt01-fwd.dcd`
   * `python xyPlot.py -i rmsd_byChain_everystep.dat -d $' \t ' -c '1,2,3,4,5,6' -m 100 -x 'time (ns)' -y 'RMSD (A)' -l "backbone;S1;S2;S3;S4;2GBI" -o rmsd_byChain_everystep.png`
   

-----------------------------------------------------------------------------------------------------------------

### `2_ligClose`
* *Description*:
   * `closeWat.tcl`
   * `lig2protDist.py`
   * `ligCloseWat.tcl`
   * `ligHbonds.tcl`
* *Usage*:
   * `vmdt -e /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/2_ligClose/closeWat.tcl -args psf dcd sel`
* *Examples*:
   * `vmdt -e /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/2_ligClose/closeWat.tcl -args ../../*.psf ../../npt01.dcd GBI2`
   * `vmdt -e /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/2_ligClose/ligHbonds.tcl -args ../../*.psf ../../npt01.dcd GBI1`


-----------------------------------------------------------------------------------------------------------------

### `3_nonCovInts`

* *Description*:
   * Analyze contact interactions between Hv1/2GBI system.
   * `getEdges.v5.tcl`
   * `getNodes.v5.tcl`
   * `nodeContacts.py`
   * `node.slurm`
* *Notes*:
   * See README inside directory for usage.
   * Original documentation: `/beegfs/DATA/mobley/limvt/hv1/04_fep/compare_r208k/nodeAnalysis/README.md`
   * See subgroup data from 2017-09-14 for R208K.i
   * See past applications:
      * `/beegfs/DATA/mobley/limvt/hv1/04_fep/compare_r208k/nodeAnalysis/`
      * `/beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/a1_D112E/1_E112/analysis/3_nonCovInts/`
      * `/beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/b1_V178A/1_A178/analysis/3_nonCovInts/` 

-----------------------------------------------------------------------------------------------------------------

### `4_dists`

* *Description*: 
   * `calcAng.tcl`
   * `calcDist.tcl` - Calculate the distance between one selection of interest and a number of other selections of interest.
* *Parameters*: 
   * PSF file
   * DCD file
   * main selection of interest, 
   * other selection to measure distance from. 
* *Returns*: 
   * Text file with info on input parameters, and columns of (1) frame number, (2) distance1, (3) distance2, etc.
   * Output file is named measDist.dat so rename if rerunning the script on new input data.
* *Usage*: 
   * `vmd -dispdev none -e calcDist.tcl -args input.psf input.dcd atom-selection1 atom-selection2 [atom-selection3] [atom-selection4]`   
* *Examples*:
   * `vmdt -e calcDist.tcl -args hHv1_open_wGBI.psf npt10.dcd protein,and,resid,211,and,name,CZ protein,and,resid,150,and,name,CG`
* *Notes*: 
   * If there are spaces in the selection, separate them with COMMAS, else will be parsed as multiple selections.
     (Quotation marks, square brackets, and curly brackets don't help.)

-----------------------------------------------------------------------------------------------------------------
