
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

04. In VMD, view name.fep file. "vmd *fep -e /beegfs/DATA/mobley/limvt/hv1/04_fep/viewMutation.tcl -args NUM GBI1"
     (-) NUM is the protein residue number.
     (a) check cell alignment with axes. see code snippets below.
     (b) Red atoms are disappearing and blue atoms are appearing (set color Beta).
     (c) Also view in Licorice "not resname GBI1 and within 3 of resname GBI1".

05. In 00_main file, edit alchemy.inp for periodic cell
     (a) cellBasisVectors can come from pdb of trajectory (not the pdb from psfgen)
     (b) ....  review all lines.

06. Edit simulation settings in slurm script.
     (a) Is the temperature correct? Is the CGenFF version correct? Is the ligand psf chg file correct?
     (b) Test a few jobs on a single node (make sure files get read in correctly, etc.)
     (c) Make sure hybrid slurm script has exclude for c-4-3 and c-4-20
     (c) Submit using hybrid slurm script on 8 nodes using "sbatch -N 8 file.slurm"

## Restarting from checkpoint

01. Check that last time step in log file matches xsc and fepout. 
    If not, modify stepline variable in checkpoint slurm script.

02. Review checkpoint slurm script. Change array numbers.

03. Submit script INSIDE the FEP_*/ directory containing lambda windows.

04. Combining windows
     (-) http://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2006-2007/3594.html
         Collect dE values from the two half-windows and combine using exponential formula
     (a) From simulation1: 
           A = value of restarted step from log file (should be consistent with xsc)
           B = last step reported in fepout file
         In lambda dir, make a copy of run1.fepout.
         Trim end lines such that last line in fepout is A, not B.
     (b) Copy all lines of run2.fepout into the copy-run1.fepout file.
     (c) Verify that lines align with a singly successfully completed run
     (d) Copy and rename as run1.fepout in the RESULTS directory.

## Data analysis

01. Python script to join windows and apply BAR.
     (a) mpython FEP_BARanalysis.py -d ../../ -e 1 -t 5 -p -v > a.out

For comparing output values of different systems / mutations, 
  * find ./ -name 'a.out' -exec grep "Net dG" {} +


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


