/* vim: set filetype=markdown : */


# Documentation for analysis scripts for FEP simulations in NAMD
Last updated: Apr 10 2018

-----------------------------------------------------------------------------------------------------------------

## Contents
* `1_dG`
* `2_rmsd`
* `3_dihed`
* `4_TTcontacts`
* `5_rmsf`
    - Use `structural/analyzeDCD.tcl`
* `6_oneWin`
* `7_waters`
    - Use `structural/analyzeDCD.tcl`
* `fixChg`
* `structural`
* `view`

-----------------------------------------------------------------------------------------------------------------

Documentation template:

## `Analysis directory`

Filename(s) go here.
* *Description*:
* *Required parameters*:
* *Returns*:
* *Usage*:
* *Example*:


-----------------------------------------------------------------------------------------------------------------


## `1_dG` analysis script

`/beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/1_dG/bar4fep.py`
* *Description*:
   * Compile .fepout files from each window of FEP calculations, then compute final free energy using BAR.
   * BAR = Bennett acceptance ratio
* *Required parameters*:
   * Directory location of where `FEP_F` and `FEP_R` are located.
   * Equil time (nanoseconds, integer)
   * Total simulation time (nanoseconds, integer)
* *Returns*:
   * `.fepout` files for each of fwd and rev results
* *Example*:
   * Execute in the main job's subdirectory of analysis results.   
   * `python /beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/1_dG/bar4fep.py -d ../../ -v -p -e 1 -t 5 --decomp > a.out`


-----------------------------------------------------------------------------------------------------------------


## `2_rmsd` analysis scripts

`/beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/2_rmsd/1_rmsdFEPwins_byChain.tcl`
* *Description*: 
   * Calculate RMSD of Hv1 protein, segregated by S1-S4 segments, over the 40 windows of FEP mutation.
* *Required parameters*:
   * Directory of FEP mutation
   * Which ligand tautomer was used
* *Returns*:
   * 
* *Usage*:
   * `vmdt -e file.tcl -args homedir withGBI`
   * homedir: full path ending with pose and mutation, e.g. /path/to/19415-10/a1_F150A
   * withGBI: 0 noGBI, 1 taut1, 2 taut2
* *Example*:  
   * `vmdt -e 1_rmsdFEPwins_byChain.tcl -args /beegfs/DATA/mobley/limvt/hv1/04_fep/0_tautomer/d4d_R208K-restr 0`
* *Notes*:
   * Should probably request an interactive job beforehand.
   * May need to change script inside if want to read both FWD and REV or just one direction only.  
     Don't leave curly brace even while commented.


-----------------------------------------------------------------------------------------------------------------


## `4_TTcontacts` analysis scripts
`/beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/4_TTcontacts/analysisInp.tcl`
`/beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/4_TTcontacts/contactAnalysis.tcl`
*From*: J. Alfredo Freites

* *Description*:
   * Evaluate contact interactions within Hv1/2GBI system.
   * Try to explain on a physical basis why binding may be preferred or not.
   * e.g., Binding seems to be less favored with this mutation because these interactions are lost.
   * [1] `analysisInp.tcl` filters contacts to those within vicinity of specified selection (e.g., 2GBI)
   * [2] `contactAnalysis.tcl` translates indices to atom names and residues and generates output files
* *Returns*:
   * [1] `analysisInp.tcl` returns a .dat file for each selection in `mySelections` with 
* *Usage*:
* *Example*:
* *Notes*:
   * Deprecated in favor of nodes/edges analysis in postFEP analysis.


-----------------------------------------------------------------------------------------------------------------


## `6_oneWin` analysis scripts

`/beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/6_oneWin/extractWork.tcl`
* *Description*: Process fepout file and associated NAMD trajectory to segregate based on cutoff of dE value of fepout.
* *Required parameters*: Input fepout file, associated PSF and DCD, filename identifier.
   * For a filename identifier of "R01", the outputs would be named `timesteps_work_R01.dat`, `alchemy_R01_bin1.fepout`, `traj_R01_bin1.dcd` 
* *Returns*:
   * dat file with list of timesteps separated into bin1 and bin2
   * bin1 fepout having snapshots where fepout dE is less than or equal to the cutoff
   * bin2 fepout having snapshots where fepout dE is less than or equal to the cutoff
   * bin1 DCD having snapshots where fepout dE is less than or equal to the cutoff
   * bin2 DCD having snapshots where fepout dE is greater than the cutoff
* *Usage*:
   * `vmdt -e extractWork.tcl -args in.fepout in.psf in.dcd outputBase`
* *Example*:
   * `vmdt -e extractWork.tcl -args ../../FEP_R/lambda_01/alchemy01.fepout ../../00_main/18629-19_R211S.psf ../../FEP_R/lambda_01/alchemy01.dcd R01`


`/beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/6_oneWin/plotMeasVsWork.py`
* *Description*:
* *Required parameters*:
* *Returns*:
* *Usage*:
* *Example*:


-----------------------------------------------------------------------------------------------------------------


## View scripts

`/beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/view/openAllWindows.tcl`
* *Description*: Open all lambda windows in FEP simulation for viewing/analyzing. 
* *Required parameters*: Location of FEP simulations, taut code (1=taut1, 2=taut2, 0=noGBI).
* *Returns*: Open VMD session with all windows loaded.
* *Usage*: `vmd -e file.tcl -args homedir withGBI`
* *Example*:
   * `vmd -e /beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/view/openAllWindows.tcl -args /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/d4d_R208K-restr 1`
* *Note*: 
   * Doesn't work on `mf_m-c1.9` or `mf_m-c2.2` partitions of interactive jobs
