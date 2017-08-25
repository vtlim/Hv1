#!/usr/bin/env python

import numpy as np
from pymbar import BAR
from pymbar import timeseries

## fwd work = work required to remove restraints

def calcWork(restrFile, unrestrFile, remove, restrSkip=0, unrestrSkip=0):
    """
    Parameters
    ----------
    restrFile   str
    unrestrFile str
    remove:     Bool
    restrSkip   int
    unrestrSkip int

    Returns
    -------
    workArray_ss

    """
    restrEne = np.loadtxt(restrFile, skiprows=restrSkip)
    unrestrEne = np.loadtxt(unrestrFile, skiprows=unrestrSkip)
    if fwd:
        work = np.subtract(unrestrEne,restrEne)
    else:
        work = np.subtract(restrEne,unrestrEne)
    return work

def subsample(workArray, label=""):
    g = timeseries.statisticalInefficiency(workArray)
    indices = timeseries.subsampleCorrelatedData(workArray, g)
    workArray_ss = workArray[indices]
    print("%s correlation time is %f" % (label, g))
    return workArray_ss


### biased energies
### unbiased energies

### ARG FWD WORK
arg_2c_yes = np.loadtxt("2a_R208-restr/PE-restr.dat", skiprows=1000)
arg_2c_not = np.loadtxt("2a_R208-restr/PE-unrestr.dat", skiprows=990)
arg_2c_fwd = np.subtract(arg_2c_yes, arg_2c_not)
#print arg_2c_yes, len(arg_2c_yes), len(arg_2c_yes[0])
#print arg_2c_not, len(arg_2c_not)
print arg_2c_fwd

### ARG REV WORK
arg_2d_yes = np.loadtxt("2b_R208-unrestr/PE-restr.dat", skiprows=990)
arg_2d_not = np.loadtxt("2b_R208-unrestr/PE-unrestr.dat", skiprows=1000)
arg_2d_rev = np.subtract(arg_2d_yes, arg_2d_not)
#print arg_2d_yes, len(arg_2d_yes), len(arg_2d_yes[0])
#print arg_2d_not, len(arg_2d_not)
print arg_2d_rev[:,1]

### BAR
lys_fwd = calcWork("2c_K208-restr/PE-restr.dat","2c_K208-restr/PE-unrestr.dat",True,2001,991)
lys_rev = calcWork("2d_K208-unrestr/PE-restr.dat","2d_K208-unrestr/PE-unrestr.dat",False,3491,4501)
lys_fwd_ss = subsample(lys_fwd[:,1])
lys_rev_ss = subsample(lys_rev[:,1])
dG, stderr = BAR(lys_fwd_ss,lys_rev_ss)
print dG, stderr


arg_fwd = calcWork("2a_R208-restr/PE-restr.dat","2a_R208-restr/PE-unrestr.dat",True,1000,990)
arg_rev = calcWork("2b_R208-unrestr/PE-restr.dat","2b_R208-unrestr/PE-unrestr.dat",False,990,1000)
arg_fwd_ss = subsample(arg_fwd[:,1])
arg_rev_ss = subsample(arg_rev[:,1])
dG, stderr = BAR(arg_fwd_ss,arg_rev_ss)
print dG, stderr
