

 ## Objective
 See how Hv1 / 2GBI behave after alchemical transformation.
 Version: Jun 07 2017

 ## Instructions
 1. Generate pdb file from last frame of fwd transformation.
    * `vmdt -e /work/cluster/limvt/analysis/writePDB.tcl -args /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/b1_V178A/00_main/18629-19_V178A.psf /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/b1_V178A/FEP_F/lambda_40/alchemy40.dcd 249 taut1_18629-19_V2A.pdb 1 1`
 
 2. Process PDB file to remove initial residue and rename final residue.
    * `mpython /work/cluster/limvt/analysis/editFEPpdb.py taut1_18629-19_V2A.pdb taut1_18629-19_V2A_fromPY.pdb V2A ALA`
 
 3. Generate PSF file and new PDB (with fixed atom numbers).
    * `vmdt taut1_18629-19_V2A_fromPY.pdb -e mem_setup.tcl > psfgen.out`
 
 4. Edit NAMD input file and run.
