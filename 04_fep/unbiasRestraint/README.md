/* vim: set filetype=markdown : */



# Description of (un)biased energies and (fwd/rev) work values
Combine energies to add or remove restraint to correct relative free energies obtained from alchemical free energy perturbation (Hv1, R208K) in NAMD.

![thermodynamic_cycle](https://cdn.pbrd.co/images/GHcQqIZ.png)


### Description of components
| dG | w/GBI | forward work | reverse work |
|---| :-: | ----------------- | --------------------- |
| 6 | no  | [1] ARG, restrained | [2] ARG, unrestrained |
| 7 | no  | [3] LYS, restrained | [4] LYS, unrestrained |
| 8 | yes | [5] ARG, unrestrained | [6] ARG, restrained |
| 9 | yes | [7] LYS, unrestrained | [8] LYS, restrained |


"You'd be calculating, say, the free energy for removing the restraints, so you would take your restrained simulation and re-evaluate the energies of that simulation without restraints (or, equivalently, just take the energies of all snapshots and subtract off the energies of the restraints at each snapshot) and the difference between those values would be your forward work. Then you'd take your unrestrained simulation and compute the energy of INTRODUCING restraints at each snapshot, and get the reverse work from that." -DLM, 07-28-2017

### Checklist
Table of index of work, dcd number of frames, colvars traj type, Greenplanet directory, and name of the .traj file.  

X = done  
O = not done  
<< = running  

| ID | dcd  | traj | dir   | trajFile |
|----| :-:  | :--: | :---- | :------- |
|  1 | 2501 | arg  | d8_2a | npt01-fwd.colvars.traj
|  2 | 2501 | arg  | d8_2b | arg-un.colvars.traj
|  3 | 2501 | lys  | d8_2c | npt01-fwd.colvars.traj
|  4 | 5001 | lys  | d8_2d | lys-un.colvars.traj
|  5 | 2500 | arg  | d7_2b | arg-un.colvars.traj
|  6 | 2500 | arg  | d7_2a | npt01.colvars.traj
|  7 | 5001 | lys  | d7_2d | lys-un.colvars.traj
|  8 | 2501 | lys  | d7_2c | npt01-fwd.colvars.traj
| 5b | 5000 | arg  | d2_2b | arg-un.colvars.traj
| 6b | 2500 | arg  | d2_2a | npt01-restr.colvars.traj
| 7b | 5001 | lys  | d2_2d | lys-un.colvars.traj
| 8b | 2843 | lys  | d2_2c | npt01-fwd.colvars.traj

### Greenplanet directories
```
[1] /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/d8_R208K-noGBI-restr/2a_R208-restr
[2] /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/d8_R208K-noGBI-restr/2b_R208-unrestr
[3] /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/d8_R208K-noGBI-restr/2c_K208-restr
[4] /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/d8_R208K-noGBI-restr/2d_K208-unrestr

taut1
[5] /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/d7_R208K-restr/2b_R208-unrestr
[6] /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/d7_R208K-restr/2a_R208-restr
[7] /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/d7_R208K-restr/2d_K208-unrestr
[8] /beegfs/DATA/mobley/limvt/hv1/04_fep/1_tautomer/18629-19/d7_R208K-restr/2c_K208-restr

taut2
[5b] /beegfs/DATA/mobley/limvt/hv1/04_fep/2_tautomer/19415-13-c1/d2_R208K-restr/2b_R208-unrestr
[6b] /beegfs/DATA/mobley/limvt/hv1/04_fep/2_tautomer/19415-13-c1/d2_R208K-restr/2a_R208-restr
[7b] /beegfs/DATA/mobley/limvt/hv1/04_fep/2_tautomer/19415-13-c1/d2_R208K-restr/2d_K208-unrestr
[8b] /beegfs/DATA/mobley/limvt/hv1/04_fep/2_tautomer/19415-13-c1/d2_R208K-restr/2c_K208-restr
```

## Procedure

1. Run all the eight MD simulations representing the mutated residue:
   * with and without the ligand
   * before and after transformation
   * with and without restraints   

   See below for example procedure for running the MD simulations.

2. For the 4 simulations that did not have restraints turned on, run the `coorfile` command on the output trajectory to get the \*.traj file with the collective variable information.
   * `cp equil01-cas.inp coorfile.inp`
   * Edit the .inp to (1) change output name, (2) include something like the following:
        ```
        ### Sidechain restraint for residue 208
        colvars          on                    
        colvarsConfig    colvars.tcl           
        
        #minimize 1000                                                   
        #run 2500000 ;# 5ns                                              
                                                                         
                                                                         
        ### Perform energy analysis without colvars on same trajectory.  
                                                                         
        set ts 1000   ; # first frame saved was frame 1000               
        coorfile open dcd XXXXXX.dcd                                      
                                                                         
        # Read all frames until nonzero is returned.                     
        while { ![coorfile read] } {                                     
          firstTimestep $ts                                              
          # Compute energies and forces, but don't try to move the atoms.
          run 0                                                          
          incr ts 1000                                                   
        }                                                                
        coorfile close                                                   
        
        ```
   * Copy over the correct colvars file (lys or arg). Verify atom numbers.
   * Run the trajectory analysis.

3. For each trajectory, extract the unaltered energies of each snapshot into a .dat file.
    * load VMD
    * `namdplot POTENTIAL vs TS out01.log` and save manually as ASCII matrix



### Sample procedure for seeding MD simulations
Source: 

