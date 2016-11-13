#!/bin/bash

#SBATCH --job-name="F150A_1504"
#SBATCH --array=1-40
#SBATCH --partition=mf_nes2.8
#SBATCH --mem=8gb
#SBATCH --time=72:00:00
#SBATCH --export=ALL


#SBATCH --ntasks-per-node=2
#SBATCH --ntasks-per-core=2
#SBATCH --threads-per-core=2

#SBATCH --cpus-per-task=8
#SBATCH --exclude="c-3-[341-344], c-6-[231-234], c-8-[251-254]"

# SBATCH --cpus-per-task=12
# SBATCH --exclude="c-4-[3-20], c-5-[9-42], c-6-[27-30], c-7-[7-24], c-8-[16-24,27,28,31,32,35,36,39,40]"

#------------SIMULATION SETTINGS -----------#
#binfile="npt29"
pdbfile="15183_04-F150A.pdb"
psffile="15183_04-F150A.psf"
fepfile="15183_04-F150A.fep"

fepsource="fep.tcl"
protpar="par_all36_prot.prm"
lippar="par_all36_lipid.prm"
cgfpar="par_all36_cgenff.prm"
watpar="toppar_water_ions.jaf2.str"
gbipar="gbi_final.str"
#restrfile="extrabonds.txt"

dLambda="0.025" # 40 Lambda Windows
temp="300"
nMinSteps=1000
eqSteps=500000  # 1 ns
nSteps=2500000  # 5 ns
#------------------------------------------#

module load namd/2.11-MPISMP
module list

cd $SLURM_SUBMIT_DIR

echo "Working Directory:" pwd
echo 'Array ID number:' $SLURM_ARRAY_TASK_ID

lambda=$SLURM_ARRAY_TASK_ID
if (( $(echo "$dLambda > 0" | bc -l) )); then
   wdir="FEP_F"
   mkdir $wdir; cd $wdir
   
   lambdalist=()
   for L in $(seq 0 $dLambda 1); do
      lambdalist+=("$L")
   done
   if [ $lambda -eq 1 ]; then
      #FEPcmd="runFEP 0 ${lambdalist[$lambda]}   $dLambda  $nSteps"
      FEPcmd="runFEPmin     0   ${lambdalist[$lambda]}   $dLambda  $nSteps  $nMinSteps $temp"
   else
      #FEPcmd="runFEP 0 ${lambdalist[$lambda]}   $dLambda  $nSteps"
      FEPcmd="runFEPmin     ${lambdalist[$lambda-1]}   ${lambdalist[$lambda]}   $dLambda   $nSteps  $nMinSteps $temp"
   fi
else
   wdir="FEP_R"
   mkdir $wdir; cd $wdir
   lambdalist=()
   for L in $(seq 1 $dLambda 0); do
      lambdalist+=("$L")
   done
   if [ $lambda -eq 1 ]; then
      FEPcmd="runFEPmin     1   ${lambdalist[$lambda]}   $dLambda  $nSteps  $nMinSteps $temp" 
   else
      FEPcmd="runFEPmin     ${lambdalist[$lambda-1]}   ${lambdalist[$lambda]}   $dLambda   $nSteps  $nMinSteps $temp" 
   fi
fi

if [ $lambda -lt 10 ]; then
   lambda_dir=`printf "lambda_0%s" $lambda`
   fname="alchemy0$lambda"
else
   lambda_dir=`printf "lambda_%s" $lambda`
   fname="alchemy$lambda"
fi

if [ ! -d "$lambda_dir" ]; then
   mkdir $lambda_dir
fi

slurm_startjob(){
cd $lambda_dir
cp ../../00_main/alchemy.inp ./$fname.inp
#Specify inputs at top
#sed -i '1iset binfile            ../../00_main/'$binfile  $fname.inp
sed -i '2iset outfile            '$fname    $fname.inp
sed -i '3istructure              ../../00_main/'$psffile  $fname.inp
sed -i '4icoordinates            ../../00_main/'$pdbfile  $fname.inp
sed -i '5iparaTypeCharmm         on' $fname.inp
sed -i '6iparameters             /data12/cmf/limvt/toppar/'$lippar  $fname.inp
sed -i '7iparameters             /data12/cmf/limvt/toppar/'$protpar   $fname.inp
sed -i '8iparameters             /data12/cmf/limvt/toppar/'$cgfpar   $fname.inp
sed -i '9iparameters             /data12/cmf/limvt/toppar/'$watpar   $fname.inp
sed -i '10iparameters             /data12/cmf/limvt/toppar/'$gbipar   $fname.inp
sed -i "11iset temp               ${temp}\n"     $fname.inp
   
#Specify FEP Parameters at end
echo "source                   ../../00_main/$fepsource" >> $fname.inp
echo "alchFile                 ../../00_main/$fepfile"   >> $fname.inp
#echo " extraBonds               on"         >> $fname.inp
#echo " extraBondsFile           ../../00_main/$restrfile" >> $fname.inp
echo "alchEquilSteps           $eqSteps"   >> $fname.inp
echo "set nMinSteps            $nMinSteps" >> $fname.inp
echo "$FEPcmd"                             >> $fname.inp
echo "#          start stop dLambda nSteps nMinSteps temp" >> $fname.inp

cpt=$SLURM_CPUS_PER_TASK
#mpiname -a


srun --propagate=STACK --cpu_bind=v,sockets --cpus-per-task=$cpt \
     namd2 ++ppn $((cpt-1)) +isomalloc_sync $fname.inp > $fname.log

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
