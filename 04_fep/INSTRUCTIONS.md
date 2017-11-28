/* vim: set filetype=markdown : */

## System setup and job submission

01. Get pdb from last frame of trajectory. In VMD, make sure that
     (a) cell is aligned with axes 'pbc box -center origin'.
     (b) membrane is parallel to xy plane, perpendicular to z axis.

02. Generate hybrid PSF:
     (a) Move psf and pdb to 01_setup directory.
     (b) In vim, edit mem_setup.tcl for the residue mutation, hybrid topology file, etc. Then run.
          * Comment GBI-related lines in the script to generate the psf (mem_setup.tcl)

03. Edit .fep file:
     (a) In 00_main dir, create copy of name.pdb to name.fep.
     (b) In vim, search for "F2A" or whatever mutation it is.
     (c) Edit next to last column - the beta values for the changing residue.
         Just because F2A does NOT mean the beta value should be changed.
         Use 3rd column as a guide, atomtypes ending with A disappear (-1.00 value), 
            atomtypes ending with B are appearing (1.00 value).
         Hint: it might help to generate a temporary empty line before & after the -1 and 1 lines. 

04. In VMD, view name.fep file. `vmd *fep -e /beegfs/DATA/mobley/limvt/hv1/04_fep/viewMutation.tcl -args NUM GBI1`
     (-) NUM is the protein residue number.
     (a) check cell alignment with axes. see code snippets below.
     (b) Red atoms are disappearing and blue atoms are appearing (set color Beta).
     (c) Also view in Licorice "not resname GBI1 and within 3 of resname GBI1".

05. In `00_main` file, edit alchemy.inp for periodic cell
     (a) cellBasisVectors can come from pdb of trajectory (not the pdb from psfgen)
     (b) ....  review all lines.

06. Edit simulation settings in slurm script.
     (a) Is the temperature correct? Is the CGenFF version correct? Is the ligand psf chg file correct?
     (b) Test a few jobs on a single node (make sure files get read in correctly, etc.)
     (c) Make sure hybrid slurm script has exclude for c-4-3 and c-4-20
     (c) Submit using hybrid slurm script on 8 nodes using "sbatch -N 8 file.slurm"

## Restarting from checkpoint

01. Check that last time step in log file matches xsc. 
    If not, modify stepline variable in checkpoint slurm script.

02. Review checkpoint slurm script. Change array numbers.

03. Submit script INSIDE the `FEP_*/` directory containing lambda windows.

04. Combining windows
     (-) http://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2006-2007/3594.html
         Collect dE values from the two half-windows and combine using exponential formula
     (a) From simulation1: 
           A = value of restarted step from log file (should be consistent with xsc)
           B = last step reported in fepout file
         In lambda dir, make a copy of run1.fepout.
         Trim end lines such that last line in fepout is A, not B.
     (b) Copy all lines of run2.fepout into the copy-run1.fepout file.
     (c) Verify same number of lines as a singly successfully completed run
     (d) Copy and rename as run1.fepout in the RESULTS directory.

## File analysis

Getting all the dG files from NAMD output
  * `grep 'Free energy' *fepout | awk '{print $12}'`

Checking out work in one particular window:
  * `xb *fepout 2 7`
  * `xmgrace -block file1.fepout -bxy 2:7 -block file2.fepout -bxy 2:7` 

For comparing output values of different systems / mutations, 
  * find ./ -name 'a.out' -exec grep "Net dG" {} +


## Restrained FEP

01. Determine appropriate restraint case to use.

02. Run FEP jobs with restraint.

03. Run restrained AND unrestrained jobs at start and final end points.   
    Also see: `/beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/README.md`  
    Subdirectories example:  
    * `2a_R208-restr`
    * `2b_R208-unrestr`
    * `2c_K208-restr`
    * `2d_K208-unrestr`

    (a) Generate pdb file from last frame of fwd transformation.  
        `vmdt -e /work/cluster/limvt/analysis/writePDB.tcl -args ../00_main/18629-19_R208K.psf ../FEP_F/lambda_40/alchemy40.dcd 249 taut1_18629-19_R2K.pdb 1 1`  
    (b) Process PDB file to remove initial residue and rename final residue.
        `python /work/cluster/limvt/analysis/editFEPpdb.py taut1_18629-19_R2K.pdb taut1_18629-19_R2K_fromPY.pdb R2K LYS`  
    (c) Generate PSF file and new PDB (with fixed atom numbers).  
        `vmdt taut1_18629-19_R2K_fromPY.pdb -e mem_setup1.tcl > psfgen.out`  
    (d) Edit NAMD input file to include RESTRAINT and run.
    (e) View...
        - After loading: `source /beegfs/DATA/mobley/limvt/hv1/04_fep/compare_r208k/zhbonds.vmd`  
        - Or this, but may need to change the vmd script to read in as pdb and not dcd:  
        `vmd -e /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/d2_R208K-noGBI/r208-general.vmd -args inpsf inpdb 1`  

04. Tcl script to calculate distances related to restraint, whether on or not.

05. Unbias endpoints using restraint definition to calculate restraint energy.

06. Correct for free energies of FEP with restraint energy.


## Miscellaneous

Job management code snippets
  * `for k in lambda_*; do echo $k; tail -n 1000 $k/*log | grep "WRITING EXTENDED SYSTEM" | tail -1; done`

VMD Code Snippets
  * protein and resid 182 and beta<=1 and beta>=0
    protein and resid 182 and beta>=-1 and beta<=0
  * measure center [atomselect top all] weight mass
  * set minmax [measure minmax [atomselect top all]]
  * molinfo top get {a b c}
  *   in graphical representations: abs(x)>48 or abs(y)>48
  * pbc box -center origin


Ligand Removal
  * vmdt -e file.pdb
  * set sel [atomselect top "not resname GBI1"]
  * $sel writepdb 15183_04_prot.pdb


Additional resources:
  * https://github.com/limn1/DTT
  * NAMD FEP tutorial PDF and website


