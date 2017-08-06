#!/usr/bin/env python

# Purpose:
# Usage:
# References:
# WARNINGS:
#  - not directly applicable to gromacs results! aside from file formatting
#    most importance difference is units. force constant in namd for angles
#    is kcal/mol/deg^2, but in gromacs is kJ/mol/rad^2


import os, glob, re
import numpy as np
from pymbar import timeseries
import pymbar


def getWindow(filename, tstart=0, tstop=0):
    """
    Read window .traj file, compute correlation times, subsample data.

    Parameters
    ----------
    filename: string name of the file to process.
       For *.traj file, assumes all lines are data (e.g. no comment lines).
    tstart: integer nanosecond start time
    tstop: integer nanosecond stop time

    Returns
    -------
    counts: int, number of entries for this particular window
    winZ: numpy list containing SUBSAMPLED data for this window from tstart to tstop

    """
    # Open file.
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()

    ### TODO: convert tstart (ns) to step number
    ### TODO: get only get data from step number and later

    # Get data from file.
    n = 0
    winZ = np.zeros([len(lines)])
    for line in lines:
        if line[0] != '#' and line[0] != '@':
            tokens = line.split()
            chi = float(tokens[1])
            # wrap ANGLE to be within [-180,+180)
            while(chi < -180.0):
                chi += 360.0
            while(chi >= +180.0):
                chi -= 360.0
            winZ[n] = chi
            n += 1

    # Compute correlation times for z (actual spring center position) timeseries.
    # use cos/sin for timeseries analysis because...?
    chi_radians = winZ[0:n]/(180.0/np.pi)
    g_cos = timeseries.statisticalInefficiency(np.cos(chi_radians))
    g_sin = timeseries.statisticalInefficiency(np.sin(chi_radians))
    print "g_cos = %.1f | g_sin = %.1f" % (g_cos, g_sin)
    g_k = max(g_cos, g_sin)
    print "Correlation time for %s =  %10.3f" % (filename,g_k)
    indices = timeseries.subsampleCorrelatedData(chi_radians, g=g_k) 

    # Subsample data.
    shortlen = len(indices)
    shortZ = winZ[indices]
    return shortlen, shortZ





# ------------------------- Parse Command Line Inputs ----------------------- #
if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option('-d','--dir',
            help = "Directory containing all colvar subdirectories.",
            type = "string",
            dest = 'hdir')

    parser.add_option('-l','--lower',
            help = "Lower bound of collective variable.",
            type = "float",
            dest = 'lower')

    parser.add_option('-u','--upper',
            help = "Upper bound of collective variable.",
            type = "float",
            dest = 'upper')

    parser.add_option('-s','--skip',
            help = "Skip incremenent of collective variable.",
            type = "float",
            dest = 'skip')
    # TODO if not evenly spaced, can't use skip or colrange

    (opt, args) = parser.parse_args()



    os.chdir(opt.hdir)
    z_min = opt.lower  # min independent variable
    z_max = opt.upper # max independent variable
    colrange = range(int(z_min), int(z_max)+int(opt.skip), int(opt.skip))
    colrange = range(0, 370, 10) # vtl fuckup

    colvarsFile = "colvars.tcl"
    trajFile = "npt01.colvars.traj"
#    trajFile = "npt01_1-5ns.traj"
    N_max = 2510 # number of snapshots max per window
    nbins = 36 # how many 'x data points' in the pmf
    
    temp = 300.
    kB = 1.381e-23 * 6.022e23 / 1000.0 / 4.184 # Boltzmann constant in kcal/mol/numWins
    beta = 1.0 / (kB * temp) # inverse temperature of simulations (in 1/(kcal/mol))
    
    
    numWins = len(colrange)
    centers = np.zeros([numWins], np.float64)
    springs = np.zeros([numWins], np.float64)
    actuals = np.zeros([numWins,N_max], np.float64)
    winLens = np.zeros([numWins], np.int32)
    
    for i, x in enumerate(colrange):
        # get center and spring constant of each window
        colvf = '%d/%s' % (x, colvarsFile)
        with open(colvf) as f:
            for line in f: 
                if 'centers' in line: centers[i] = line.split()[1]
                if 'forceConstant' in line: springs[i] = line.split()[1]
    
    
        # process simulation data for each window
        trajf = '%d/%s' % (x, trajFile)
        with open(trajf) as f:
            ilen, winZ = getWindow(trajf)
            actuals[i][0:ilen] = winZ
            winLens[i] = ilen
   
    # vtl fuckup wrap centers [-180,+180) ================================ ***
    for i, chi in enumerate(centers):
        print chi
        if chi < -180.0:
            centers[i] += 360.0
        if chi >= +180.0:
            centers[i] -= 360.0

    N_max = np.max(winLens)
    u_kn = np.zeros([numWins,N_max], np.float64)
    u_kln = np.zeros([numWins,numWins,N_max], np.float64)
    
    # Construct torsion bins
    print "Binning data..."
    delta = (z_max - z_min) / float(nbins)
    # compute bin centers
    bin_center_i = np.zeros([nbins], np.float64)
    for i in range(nbins):
        bin_center_i[i] = z_min + delta/2 + delta * i
    # Bin data
    bin_kn = np.zeros([numWins,N_max], np.int32)
    for k in range(numWins):
        for n in range(winLens[k]):
            # Compute bin assignment.
            bin_kn[k,n] = int((actuals[k,n] - z_min) / delta)
    
    # Evaluate reduced energies in all umbrellas
    print "Evaluating reduced potential energies..."
    for k in range(numWins):
        for n in range(winLens[k]):
            # Compute minimum-image torsion deviation from umbrella center l
            dchi = actuals[k,n] - centers
            for l in range(numWins):
                if (abs(dchi[l]) > 180.0):
                    dchi[l] = 360.0 - abs(dchi[l])
    
            # Compute energy of snapshot n from simulation k in umbrella potential l
            u_kln[k,:,n] = u_kn[k,n] + beta * (springs/2.0) * dchi**2

    # Initialize MBAR.
    print "Running MBAR..."
    mbar = pymbar.MBAR(u_kln, winLens, verbose = True, method = 'adaptive')
    
    # Compute PMF in unbiased potential (in units of kT).
    print bin_kn
    (f_i, df_i) = mbar.computePMF(u_kn, bin_kn, nbins)
    
    # Write out PMF
    print "PMF (in units of kT)"
    print "%8s %8s %8s" % ('bin', 'f', 'df')
    for i in range(nbins):
        print "%8.1f %8.3f %8.3f" % (bin_center_i[i], f_i[i], df_i[i])
    
    # Output, relabeling bin angles and energies in kcal/mol
    print "\nPMF (in units of kcal/mol)"
    print "%8s %8s %8s" % ('bin', 'f', 'df')
    for i in range(nbins):
        print "%8.1f %8.3f %8.3f" % (bin_center_i[i], f_i[i]*0.593, df_i[i]*0.593)
