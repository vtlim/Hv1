#!/usr/bin/env bash

# Things to check:
#   1. Psi has completed successfully.
#   2. Optimization completed successfully.
#   3. Optimization failed.

wdir='/work/cluster/limvt/hv1/05_dihedScan/2_tautomer/dihed-180'
cd $wdir

### Check all conformers with MP2 method
tail -n 1 ./*/*/output.dat > $wdir/did-psiFinish
grep "* Optimization is complete" ./*/*/output.dat > $wdir/did-optFinish-yes
grep "* Optimization has failed" ./*/*/output.dat > $wdir/did-optFinish-no

### Check only conformer 1, all methods
#tail -n 1 AlkEthOH*/1/*/*/output.dat > $wdir/conf1-did-psiFinish
#grep "* Optimization is complete" AlkEthOH*/1/*/*/output.dat > $wdir/conf1-did-optFinish-yes
#grep "Optimizer: Did not converge" AlkEthOH*/1/*/*/output.dat > $wdir/conf1-did-optFinish-no


#####
# First process 'did-psiFinish' in vim and rename file 'did-psiFinish-no'
# in vim, (:) g/beer/-1 d 3
# a result of "no lines in buffer" is good! psi at least terminated normally.

# Total number of jobs = 'did-psiFinish-no' + 'did-optFinish-yes' + 'did-optFinish-no'
# which can be checked with something like this:
# find ./ -mindepth 2 -type f -name 'output.dat' -exec tail {} \; 
