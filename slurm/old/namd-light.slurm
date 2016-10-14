#!/bin/bash

#SBATCH --job-name="ref_F150A"
#SBATCH --partition=mf_ilg2.3

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --distribution=block:cyclic
#SBATCH --ntasks-per-core=2
#SBATCH --threads-per-core=2
#SBATCH --exclude=c-16-24
#SBATCH --array=1
#SBATCH --mem=8gb
#SBATCH --time=72:00:00


#------------------------------------------#

module load namd/2.11-OpenMPI_gcc_ilg
module list

cd $SLURM_SUBMIT_DIR

echo "Working Directory:" pwd
echo 'Array ID number:' $SLURM_ARRAY_TASK_ID

lambda=$SLURM_ARRAY_TASK_ID

if [ $lambda -lt 10 ]; then
   lambda_dir=`printf "lambda_0%s" $lambda`
   fname="alchemy0$lambda"
else
   lambda_dir=`printf "lambda_%s" $lambda`
   fname="alchemy$lambda"
fi


slurm_startjob(){
cd $lambda_dir

HOSTFILE=./hosts.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

srun hostname -s > $HOSTFILE
opt="--hostfile $HOSTFILE
     --bind-to core
     --map-by core:pe=1
     --report-bindings
     -mca btl openib,sm,self"

cpuset=$(cat /proc/self/status | grep Cpus_allowed_list | awk '{print $2;}')
echo "Rank $OMPI_COMM_WORLD_RANK bound to core(s) $cpuset"
     
mpirun $opt namd2 +setcpuaffinity +isomalloc_sync $fname.inp > $fname.log

cd ..

echo "JOB DONE"
}

# Function to echo informational output
slurm_info_out(){
# Informational output
echo "=================================== SLURM JOB ==================================="
date
echo
echo "The job will be started on the following node(s):"
echo $SLURM_JOB_NODELIST
echo
echo "Slurm User:         $SLURM_JOB_USER"
echo "Run Directory:      $(pwd)"
echo "Job ID:             ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Job Name:           $SLURM_JOB_NAME"
echo "Partition:          $SLURM_JOB_PARTITION"
echo "Number of nodes:    $SLURM_JOB_NUM_NODES"
echo "Number of tasks:    $SLURM_NTASKS"
echo "Submitted From:     $SLURM_SUBMIT_HOST:$SLURM_SUBMIT_DIR"
echo "=================================== SLURM JOB ==================================="
echo
echo "--- SLURM job-script output ---"
}

copy_results(){

if [ ! -d results ]; then
   mkdir results
fi
cp $lambda_dir/$fname.fepout ./results/
}

slurm_startjob
copy_results
slurm_info_out

date
