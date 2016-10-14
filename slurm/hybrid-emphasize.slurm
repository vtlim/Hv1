#!/bin/bash

#SBATCH --job-name="F182Ap_fwd"
#SBATCH --array=1-24
#SBATCH --partition=mf_nes2.8
#SBATCH --mem=8gb
# SBATCH --time=72:00:00
#SBATCH --time=2:00:00
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
pdbfile="15183_04-F182A.pdb"
psffile="15183_04-F182A.psf"
fepfile="15183_04-F182A.fep"

fepsource="fep.tcl"
protpar="par_all36_prot.prm"
lippar="par_all36_lipid.prm"
cgfpar="par_all36_cgenff.prm"
watpar="toppar_water_ions.jaf2.str"
gbipar="gbi_final.str"
#restrfile="extrabonds.txt"

dLambda="0.05" # 19 Lambda Windows
dLastLam="0.01" # break up last window (1-dLambda)
temp="300"
nMinSteps=1000
eqSteps=500000  # 1 ns
nSteps=5500000  # 11 ns
#------------------------------------------#

### Comments to code explained using an example of 20 normal windows,
### with the last one broken up to 5, leading to 19+5=24 windows.
### This means dLambda = 0.05 and dLastLam = 0.01.
### Note that sign of dLambda and dLastLam should be the same.

module load namd/2.11-MPISMP
module list

cd $SLURM_SUBMIT_DIR

echo "Working Directory:" pwd
echo 'Array ID number:' $SLURM_ARRAY_TASK_ID

window=$SLURM_ARRAY_TASK_ID

### Get absolute value by removing any minus sign from dLambda
absDL=${dLambda//-}
### Get value of what would be the normal last window (e.g. 20.000).
normlast=$(echo "1 / $absDL" | bc -l)
### Get value of lambda at the start of the last window (e.g. 0.95).
lastlambda=$(echo "1-$absDL" | bc -l)

### For positive dLambda, set up forward calculations.
if (( $(echo "$dLambda > 0" | bc -l) )); then
   wdir="FEP_F"
   mkdir $wdir; cd $wdir

   ### Break the last window into smaller windows, by dLastLam. 
   if (( $(echo "$window >= $normlast" | bc -l) )); then
       ### Get list of lambda values for last window (e.g. 0.95 0.96 0.97 0.98 0.99 1.00)
       lastlist=()
       for K in $(seq $lastlambda $dLastLam 1); do
          lastlist+=("$K")
       done

       ### Get index of correct lambda to use in lastlist.
       ### E.g. if window is 24, lindex is 24-20+1 = 4, and lastlist[4]=0.99.
       lindex=$(echo "$window - $normlast + 1" | bc -l) 
       lindex=${lindex%.*} # remove trailing zeroes
       FEPcmd="runFEPmin     ${lastlist[$lindex-1]}   ${lastlist[$lindex]}   $dLastLam   $nSteps  $nMinSteps $temp"

   ### If window number < normlast, treat as normal.
   else
       lambdalist=()
       for L in $(seq 0 $dLambda $lastlambda); do
          lambdalist+=("$L")
       done
       FEPcmd="runFEPmin     ${lambdalist[$window-1]}   ${lambdalist[$window]}   $dLambda   $nSteps  $nMinSteps $temp"
   fi

### For negative dLambda, set up reverse calculations.
else
   wdir="FEP_R"
   mkdir $wdir; cd $wdir

   ### Number of subwindows = dLambda / dLastLambda
   numSubWin=$(echo "$dLambda / $dLastLam" | bc -l)
   numSubWin=${numSubWin%.*} # remove trailing zeroes

   ### E.g. The first five windows of rev ~ the last five windows of fwd.
   if [ $window -le $numSubWin ]; then
       ### Get list of lambda values for last window (e.g. 1.00 0.99 0.98 0.97 0.96 0.95)
       lastlist=()
       for K in $(seq 1 $dLastLam $lastlambda); do
          lastlist+=("$K")
       done
       FEPcmd="runFEPmin     ${lastlist[$window-1]}   ${lastlist[$window]}   $dLastLam   $nSteps  $nMinSteps $temp"

   else
       lambdalist=()
       for L in $(seq $lastlambda $dLambda 0); do
          lambdalist+=("$L")
       done

       ### Get index of correct lambda to use in lastlist.
       ### E.g. if window is 6, lindex is 6-5-1 = 0, and lambdalist[0]=0.95.
       lindex=$(echo "$window - $numSubWin - 1" | bc -l) 
       lindex=${lindex%.*} # remove trailing zeroes
       FEPcmd="runFEPmin     ${lambdalist[$lindex]}   ${lambdalist[$lindex+1]}   $dLambda   $nSteps  $nMinSteps $temp" 
   fi
fi

### Get names for window directory and input files.
if [ $window -lt 10 ]; then
   lambda_dir=`printf "lambda_0%s" $window`
   fname="alchemy0$window"
else
   lambda_dir=`printf "lambda_%s" $window`
   fname="alchemy$window"
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
