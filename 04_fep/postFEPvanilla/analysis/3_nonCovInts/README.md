
# Quantify contacts within protein-ligand system
Loosely based on node definitions by [Benson and Daggett](https://dx.doi.org/10.1142%2FS0219720012500084).
Tcl scripts written by Eric Wong and Alfredo Freites.
Documentation last updated: *Mar 23 2018*.

## Script sources
Copy original from HPC: `crsync sshhpc:/pub/limvt/hv1/02_configs/02_analysis/07_nodeContacts/*tcl /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/3_nonCovInts/`
Copy template for here: `crsync /beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/3_nonCovInts/*tcl .`

## Notes
* Directory names in the Tcl scripts should end with /
* Running Tcl scripts are pretty fast but the python script is resource-intensive

* submitting as a slurm job takes 3-4 days for 500 frames for 4 simulations
* if requesting interactive job, make sure you request appropriate time (if x hours for 1 set of 500 then 3x hours for 3 sets of 500)
* maybe turn off printing when submitting job
* maybe submit job with copy-local option

## Procedure
`/beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/3_nonCovInts/README.md`

0. Make subdirectories in the post-FEP-trajectory folder.
    * `cd /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/f1_R211S/1_S211`
    * *`mkdir -p analysis/1_rmsd analysis/2_ligClose analysis/3_nonCovInts`*
    * Copy the .tcl scripts and this procedure README to the nonCovInts directory
      (don't make a typo and forget the destination directory else you overwrite the node script!)

0. Request an interactive job (not the qrsh-mc19).

1. Acquire node network file with `getNodes.v5.tcl`. Use with *`vmdt -e getNodes.v5.tcl`*
    * Variables to change before using: [1] input/output directories, [2] PSF filename, [3] ligand residue
    * DON'T FORGET THE LIGAND RESIDUE NAME
    * For case with noGBI, set ligand residue "" which will throw error below but output is fine
      `can't read "ndict()": no such element in array`
    * Output is a bonds file and a nodes file based off the protein structure.
    * Format of nodes file: nodeIndex, resname, resid, sidechain/backbone/water/ligand, node type (positive, nonpolar, etc), nodename(?), numAtoms, molIndex

2. Create a wrapped summary trajectory.
    * Parameters: `vmdt -e sumtraj.tcl -args psffile outputName skip dcdfile1 lastframe1 [dcdfile2 lastframe2]`
    * Example:    *`vmdt -e /data12/cmf/limvt/analysis/sumtraj.tcl -args ../../taut1_18629-19_s211.psf ../../wrap10.dcd 10 ../../npt01-fwd.dcd 2501`*

3. Acquire edge information with `getEdges.v5.tcl`. Use with *`vmdt -e getEdges.v5.tcl`*
   Variables to change before using:
    1. input/output directories
    2. PSF filename
    3. trajectory (list, if 2+)
    4. trajectory range (list if 2+) (leave empty "" if analyzing all frames)
    5. selection sentence for analysis (DON'T FORGET THE LIGAND RESIDUE NAME)
    6. output file name!!
    7. start frame number of output

4. Process node edge results using `/beegfs/DATA/mobley/limvt/hv1/04_fep/postFEPvanilla/analysis/3_nonCovInts/nodeContacts.py`.
    * Do this as a slurm job bc intensive (time and memory).
    * Example: *`python nodeContacts.py -n hHv1_open_wGBI.psf.nodes -e t1-18629-19-npt0910 -a 001 -b 500 -o t1ref.pickle &`*
    * Python documentation:
```
$ python nodeContacts.py -h
usage: nodeContacts.py [-h] -n N -e E -a A -b B -o O

optional arguments:
  -h, --help  show this help message and exit
  -n N        Name of the node file. Suffix should be .nodes
  -e E        Prefix of the edge file before frame number. Don't include '_'
  -a A        First frame number of edge file to read. Include leading zeroes
              if present.
  -b B        Last frame number of edge file to read.
  -o O        Name of output file with pickled data.
```

5. Analyze node edge results in iPython notebook and plot data.
