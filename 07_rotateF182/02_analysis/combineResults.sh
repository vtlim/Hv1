
# cd to the directory with angle subdirectories
# then call this script. check relative placement of
# results directory is correct. (the line with ln -s)
# TODO command line arguments

MYPWD=$(pwd)
echo $MYPWD
for i in `ls -d */`; do
  if [ -f ${i}/npt01.colvars.traj ]; then
    echo "npt01.colvars${i%/}.traj"
    #ln -s ${MYPWD}/${i}/*traj ../results5/npt01.colvars${i%/}.traj
    ln -s ${MYPWD}/${i} ../results2/${i%/}
  fi
done
