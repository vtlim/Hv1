#!/bin/bash

# SUBMIT THIS SCRIPT INSIDE THE FEP_*/ DIRECTORY CONTAINING LAMBDA WINDOWS
# Change array number to submit individual jobs.
# If you plan to change margin (or other) do so manually in *.inp file first.

#SBATCH --job-name="FEP_test8"
#SBATCH --array=8

#SBATCH --partition=mf_nes2.8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --time=72:00:00
#SBATCH --distribution=block:cyclic


#------------SIMULATION SETTINGS -----------#
lambda=$SLURM_ARRAY_TASK_ID
#------------------------------------------#

printf "Start Time:$( date )\n"
module load namd/2.11-OpenMPI_gcc

cd $SLURM_SUBMIT_DIR
echo "Hostname:" $HOSTNAME
echo "Working Directory: $( pwd )"
echo "Array ID number:" $SLURM_ARRAY_TASK_ID

# set working lambda directories
if [ $lambda -lt 10 ]; then
   lambda_dir=`printf "lambda_0%s" $lambda`
   fname="alchemy0$lambda"
else
   lambda_dir=`printf "lambda_%s" $lambda`
   fname="alchemy$lambda"
fi

# check output file
finline=$( tail -5 "$lambda_dir"/"$fname.log" | grep "Program finished" )
if [[ $finline == *"Program finished"* ]]; then
   echo "$lambda_dir Complete"
   echo "End Time:$( date )"
   exit 
else
   echo "$lambda_dir Incomplete"
   cd $lambda_dir
   pwd
   
   # Look for restart files
   err=0

   bincoor="$fname.coor"
   if [ ! -f $bincoor ]; then
      echo "Could not find bincoor file."
      (( err++ ))
   fi
   binvel="$fname.vel"
   if [ ! -f $binvel ]; then
      echo "Could not find binvel file."
      (( err++ ))
   fi
   binxsc="$fname.xsc"
   if [ ! -f $binxsc ]; then
      echo "Could not find binxsc file."
      (( err++ ))
   fi

   # Get last timestep from which to restart
   stepline=$( tail -n 1000 "$fname.log" | grep "WRITING EXTENDED SYSTEM" | tail -1 )
   if [ -z "$stepline" ]; then
      echo "Could not find last timestep."
      (( err++ ))
   fi

   # Quit if restart files are not found
   if [ $err == 4 ]; then
      echo "Could not find checkpoint files. Restart the run manually."
      exit
   fi

   # Get the last time step
   # valid if outputfreq matches restartfreq
   echo "Editing input files to restart from checkpoint..."
   echo "Found: $stepline"
   step=${stepline##* }
   if ! [[ "$step" =~ ^[0-9]+$ ]]; then
      echo "Last timestep was not an integer. Check input file."
      exit 
   fi
   echo "Last timestep: $step"

   # Adjust input file to restart FEP from checkpoint
   cp $fname.inp $fname-cpt.inp
   sed -i "s/$fname/$fname-cpt/"               $fname-cpt.inp # outfile name* b4 input
   sed -i "/set binfile/d"                     $fname-cpt.inp # remove old inputname
   sed -i "1iset binfile            $fname"    $fname-cpt.inp # set new input
   sed -i "/#firsttimestep/firsttimestep $step" $fname-cpt.inp # set firsttimestep
   sed -i "s/temperature\t/\#temperature\t/"   $fname-cpt.inp # remove temperature* \t
   sed -i "s/#bincoordinates/bincoordinates/"  $fname-cpt.inp # use binary coordinates
   sed -i "s/#binvelocities/binvelocities/"    $fname-cpt.inp # use binary velocities
   sed -i "s/#extendedSystem/extendedSystem/"  $fname-cpt.inp # use periodic cell
   sed -i "s/cellBasis/\#cellBasis/"            $fname-cpt.inp # remove cell sizes
   sed -i "s/cellOrigin/\#cellOrigin/"          $fname-cpt.inp # remove cell origin
#   sed -i "s/alchEquil/\#alchEquil/"            $fname-cpt.inp # remove FEP Equil steps =================
        
   # Evaluate remaining FEP steps
   fepline=$( tail -n 5 "$fname.inp" | grep "runFEP" )
   if [ -z "$fepline" ]; then
      echo "Could not find total steps. Check input file."
      exit 
   fi
   echo "Found: $fepline"

   # Get the original total amount of steps
   totalstep=$( echo $fepline | awk '{print $5}' | bc )
   if ! [[ "$totalstep" =~ ^[0-9]+$ ]]; then
      echo "Total step was not an integer. Check input file."
      exit 
   fi   
   echo "Total Steps: $totalstep"

   # Subtract totalstep - laststep
   let remstep="$totalstep-$step"
   if ! [[ "$remstep" =~ ^[0-9]+$ ]]; then
      echo "Remaining steps was not an integer."
      exit 
   fi

   echo "Remaning FEP steps: $remstep"
   # runFEPmin to runFEP, keep first 35 chars, include remstep 
   sed -i -e "s/Pmin/P/" -r -e "s/runFEP(.{29}).*/runFEP\1 $remstep/" $fname-cpt.inp
        
   # Submit FEP run
   # Use OpenFabrics/IB network for internode MPI communications and
   #    shared memory for intranode MPI communications.
   # export OMP_NUM_THREADS=2
   # ompi_info
   HOSTFILE=./hosts.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
   
   srun hostname -s > $HOSTFILE
   opt="--hostfile $HOSTFILE
        --bind-to core
        --map-by core:pe=1
        --report-bindings
        -mca btl openib,sm,self"
   
   cpuset=$(cat /proc/self/status | grep Cpus_allowed_list | awk '{print $2;}')
   echo "Rank $OMPI_COMM_WORLD_RANK bound to core(s) $cpuset"
        
   mpirun namd2 +isomalloc_sync $fname-cpt.inp > $fname-cpt.log
   cd ..
fi
echo "End Time:$( date )"
exit
