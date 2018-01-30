

# Corrections for free energy calculations involving change of charge

Last updated: Dec 15 2017 

The contents here are written by [Nathan Lim](https://github.com/nathanmlim/DTT/tree/master/shell/solve-pbeq) and adapted for use for Hv1 system by Victoria Lim. 
Analytical corrections implemented in accordance with [Rocklin paper](https://doi.org/10.1063/1.4826261).

Also see:

## Execution

1. Go to FEP directory.
2. Create subdirectories.
3. Navigate to first set of calculations.
4. Copy slurm script.
5. Edit slurm script parameters.
6. Submit.

    * `mkdir -p chgcorr/FEP_F/PB_01 chgcorr/FEP_F/PB_40 chgcorr/FEP_R/PB_01 chgcorr/FEP_R/PB_40`
    * `cd /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/f1_R211S/chgcorr/FEP_F/PB_01`

    * `cp /beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/fixChg/solve-pbeq-t1.sh .`
    * Common changes within an FEP dir is FEP folder, resname, window number.

    * `vi solve-pbeq-t1.sh`
    * `sb solve-pbeq-t1.sh`

There's going to be a slurm output file for every frame (250 frames * 4 PBs * 2 directions * 3 tauts).
```
mkdir outSlurm
mv slurm*out outSlurm
tar -zcvf outSlurm.tar.gz outSlurm
[to open, untested] tar -xzf --to-command=cat file-in-archive
```


## Documentation


**`solve-pbeq.sh`**
* Description:
   * Master script; built for execution on Greenplanet with SLURM job manager.
* Parameters


**`01_dcd-crd.tcl`**
* Description: 
   * Processes NAMD trajectory for subsequent reading into CHARMM.
   * Works by wrapping system, writing out PDB for each frame, then generating CRD files for each of those PDBs.
   * Called in the `prep4charmm` bash function.
* Parameters:
   1. Frame number of analysis
   2. PSF file name
   3. FEP file name (PDB format with hybrid residue)
   4. DCD file name
   5. Location of output for PDB file
   6. Location of output for CRD file
* Returns:
   1. PDB file of the ith frame in the dcd with -----?
   2. Corresponding CRD file to that PDB file.
* Check/modify before using: 
   0. ligand selection in all of the below (in addition to Hv1 selection)
   1. `sel` variable for residue number of the mutation (2 changes, resid & greater/less than sign)
   2. line under "Restore orig. resname"                (2 changes, resnumber & resname [for start or end residue])
   3. `beta` variable                                   (2 changes, resid & greater/less than sign)
   4. `atomselect` line in foreach loop                 (1 change of resid)
   5. `pbc wrap` command                                (1 change of resid)
* Notes:
   * Can view CRD file in VMD if loaded in under correct format (?). Or copy .crd file to a .cor extension.


**`02_gencharmm.inp`**
* Description: This takes in the NAMD CRD file generated in **`01_dcd-crd.tcl`** and uses it to generate a PSF file and a CRD file for CHARMM.
* Parameters: 
   1. title
   2. i
   3. sel
   4. incrd
   5. ref
   6. toppar
   7. outcrd
   8. outpsf
* Returns:
* Check/modify before using:
   1. Protonation states of acidic residues. Something like this VMD selection: `sidechain and resname ASP GLU`. Add lines to patch these if necessary.
Check the CHARMM generated CRD file to make sure it's consistent.
* Notes:
   1. 
