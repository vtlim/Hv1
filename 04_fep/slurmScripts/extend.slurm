#!/bin/bash

# SUBMIT THIS SCRIPT INSIDE THE FEP_*/ DIRECTORY CONTAINING LAMBDA WINDOWS
# Change array number to submit individual jobs.
# If you plan to change margin (or other) do so manually in *.inp file first.

#SBATCH --job-name="t1_r2s_fwd-40"
#SBATCH --array=40
#SBATCH --partition=mf_nes2.8
#SBATCH --mem=8gb
#SBATCH --time=72:00:00
#SBATCH --export=ALL
#SBATCH --ntasks-per-node=2
#SBATCH --ntasks-per-core=2
#SBATCH --threads-per-core=2
#SBATCH --cpus-per-task=8
#SBATCH --exclude="c-5-27,c-4-3,c-4-20,c-3-[341-344], c-6-[231-234], c-8-[251-254]"

# SBATCH --cpus-per-task=12
# SBATCH --exclude="c-4-[3-20], c-5-[9-42], c-6-[27-30], c-7-[7-24], c-8-[16-24,27,28,31,32,35,36,39,40]"

#------------SIMULATION SETTINGS -----------#
lambda=$SLURM_ARRAY_TASK_ID
totalstep=5000000
#------------------------------------------#

printf "Start Time:$( date )\n"
module load namd/2.11-MPISMP
module list

cd $SLURM_SUBMIT_DIR
echo "Hostname:" $HOSTNAME
echo "Working Directory: $( pwd )"
printf "\nArray ID number:" $SLURM_ARRAY_TASK_ID

# set working lambda directories
if [ $lambda -lt 10 ]; then
   lambda_dir=`printf "lambda_0%s" $lambda`
   fname="alchemy0$lambda"
else
   lambda_dir=`printf "lambda_%s" $lambda`
   fname="alchemy$lambda"
fi
cd $lambda_dir


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
#sed -i -e "s/#firsttimestep/firsttimestep $step/" $fname-cpt.inp # set firsttimestep, not allowed for fep
sed -i "s/temperature\t/\#temperature\t/"   $fname-cpt.inp # comment temperature* \t
sed -i "s/#bincoordinates/bincoordinates/"  $fname-cpt.inp # use binary coordinates
sed -i "s/#binvelocities/binvelocities/"    $fname-cpt.inp # use binary velocities
sed -i "s/#extendedSystem/extendedSystem/"  $fname-cpt.inp # use periodic cell
sed -i "s/cellBasis/\#cellBasis/"            $fname-cpt.inp # comment cell sizes
sed -i "s/cellOrigin/\#cellOrigin/"          $fname-cpt.inp # comment cell origin
sed -i "s/alchEquil/\#alchEquil/"            $fname-cpt.inp # comment FEP Equil steps
        
# Evaluate remaining FEP steps
fepline=$( tail -n 5 "$fname.inp" | grep "runFEP" )
if [ -z "$fepline" ]; then
   echo "Could not find total steps. Check input file."
   exit 
fi
echo "Found: $fepline"

# Subtract totalstep - laststep
let remstep="$totalstep-$step"
if ! [[ "$remstep" =~ ^[0-9]+$ ]]; then
   echo "Remaining steps was not an integer."
   exit 
fi

printf "Remaining FEP steps: $remstep\n\n"
# runFEPmin to runFEP, keep first 28 chars, include remstep 
sed -i -e "s/Pmin/P/" -r -e "s/runFEP(.{22}).*/runFEP\1 $remstep/" $fname-cpt.inp
sed -i "/dLambda/d" $fname-cpt.inp
echo "#          start stop dLambda nSteps" >> $fname-cpt.inp
     
# Submit FEP run
cpt=$SLURM_CPUS_PER_TASK
#mpiname -a


srun --propagate=STACK --cpu_bind=v,sockets --cpus-per-task=$cpt \
     namd2 ++ppn $((cpt-1)) +isomalloc_sync $fname-cpt.inp > $fname-cpt.log
cd ..

echo "End Time:$( date )"
exit
