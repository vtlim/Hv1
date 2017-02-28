#!/bin/sh

# Purpose: set up subdirectories and input files for umbrella sampling.
# Usage: "./file.sh min max incr"
# Example: "file.sh -130 70 10"
#
# Min/max = min/max values for "centers" in colvars file 
# incr    = incremental step going from min to max
#
# Before using, you should have: 
#  * Base namd input file
#  * Base colvars input file
#     > with general values for colvarstrajfrequency, forceConstant, etc.
#  * Base job submission script
#  * Main directory that will house all windows
#
# You can make the colvars and output files different names
#   for each window, but will make file processing a bit harder.

# ======================================

cd /pub/limvt/hv1/07_rotateF182/windows

qsub=/pub/limvt/hv1/07_rotateF182/01_setup/win-namd.sh
colvprod=/pub/limvt/hv1/07_rotateF182/01_setup/win-colv.tcl
prodinp=/pub/limvt/hv1/07_rotateF182/01_setup/win-equil.inp

# ======================================

minC=$1
maxC=$2
incr=$3

curC=$minC  # current center value for looping
while [ "$curC" -le "$maxC" ]; do
  echo $curC

  # make subdirectories
  subd=angle${curC}
  if [ ! -d "$subd" ]; then
    mkdir "$subd"
  else
    echo ">>> Directory $subd already exists, skipping."
    curC=`expr "$curC" + "$incr"`;
    continue # don't do the below if subdir exists
  fi

  # copy namd submission file to all production directories
  sed -i "/#$ -N/c #$ -N dihed${curC}" $qsub
  cp $qsub "$subd/namd-multicore.sh"

  # edit and copy production files
  sed -i "s/centers.*/centers        $curC/g" $colvprod
  cp $colvprod "$subd/colvars.tcl"

  cp $prodinp "$subd/equil1.inp"

  curC=`expr "$curC" + "$incr"`;
done

