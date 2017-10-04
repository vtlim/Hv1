

Version: Oct 14 2017
Source:  `/beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/README.md` 

## Objective

See how Hv1 / 2GBI behave after alchemical transformation.

## Instructions for generating post-FEP vanilla trajectory

1. Generate pdb file from last frame of fwd transformation.
    * `vmdt -e /work/cluster/limvt/analysis/writePDB.tcl -args /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/b1_V178A/00_main/18629-19_V178A.psf /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/b1_V178A/FEP_F/lambda_40/alchemy40.dcd 249 taut1_18629-19_V2A.pdb 1 1`
 
2. Process PDB file to remove initial residue and rename final residue.
    * `python /work/cluster/limvt/analysis/editFEPpdb.py taut1_18629-19_V2A.pdb taut1_18629-19_V2A_fromPY.pdb V2A ALA`
 
3. Generate PSF file and new PDB (with fixed atom numbers).
    * `vmdt taut1_18629-19_V2A_fromPY.pdb -e mem_setup.tcl > psfgen.out`
 
4. Edit NAMD input file and run.

## Analyses
Do not copy over scripts except for the two Tcl scripts in nonCovInts directory. 
Otherwise call scripts by name and command line arguments so all analyses are uniform.
Update this README: by uncommenting, updating filenames and directories.

**View contacts to some selection:**
```
#vmd -e /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/viewHbonds.vmd -args ../taut1_18629-19_k208restr.psf ../npt01-fwd.dcd
```

###  `1_rmsd`

```
#vmdt -e /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/1_rmsd/rmsd_trajectory.tcl -args ../../../00_main/18629-19_V178A.pdb 1 ../../taut1_18629-19_a178.psf ../../npt01.dcd
#python /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/1_rmsd/xyPlot.py -i rmsd_byChain_everystep.dat -d $' \t ' -c '1,2,3,4,5,6' -m 100 -x 'time (ns)' -y 'RMSD (A)' -l "backbone;S1;S2;S3;S4;2GBI" -o rmsd_byChain_everystep.png
```

###  `2_ligClose`

**ligand hbond contacts**
```
#vmdt -e /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/2_ligClose/ligHbonds.tcl -args ../../*.psf ../../npt01.dcd GBI1
```

**waters around ligand**
```
#vmdt -e /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/2_ligClose/closeWat.tcl -args ../../*.psf ../../npt01.dcd GBI2
```

###  `3_nonCovInts`

Ingredients:
* sumtraj.tcl
* getNodes.v5.tcl
* getEdges.v5.tcl

Instructions, with examples:

0. Request an interactive job (not the qrsh-mc19).

1. Create a wrapped summary trajectory.
    * `#vmdt -e /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/sumtraj.tcl -args d8_2d/noGBI_k208restr.psf 10 d8_2d/wrap10.dcd d8_2d/npt01-fwd.dcd 5000`
    * ``

2. Change variables then acquire node network file with `getNodes.v5.tcl`. Use with `vmdt -e getNodes.v5.tcl`  
    * input/output directories
    * PSF file variable
    * Ligand residue

   ```
   set workDir /beegfs/DATA/mobley/limvt/hv1/04_fep/2_tautomer/19415-13-c1/d2_R208K-restr/2d_K208-unrestr/analysis/3_nonCovInts/
   set dataDir /beegfs/DATA/mobley/limvt/hv1/04_fep/2_tautomer/19415-13-c1/d2_R208K-restr/2d_K208-unrestr/
   set myPSF taut2_19415-13-c1_k208restr.psf
   ```

3. Change variables then acquire edge information with `getEdges.v5.tcl`. Use with `vmdt -e getEdges.v5.tcl`
    * input/output directories
    * PSF file
    * list of trajectories
    * list of trajectory ranges
    * selection sentence for analysis
    * output file name
    * start frame number of output

4. Process node edge results using `/beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/3_nonCovInts/nodeContacts.py`.
   ```
   #python /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/3_nonCovInts/nodeContacts.py -n taut1_npt0910/hHv1_open_wGBI.psf.nodes -e taut1_npt0910/t1-18629-19-npt0910 -a 001 -b 500 -o taut1_npt0910/t1ref.pickle &
   ```
