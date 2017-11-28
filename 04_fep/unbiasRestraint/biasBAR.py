#!/usr/bin/env python

# TODO: update plot to have horizontal lines for guides
# TODO: update file reader to read multiple (3+) columns
# TODO: make equil discard time an argument
# TODO: incorporate equil discard time by not plotting that

import sys
import numpy as np
from pymbar import BAR
from pymbar import timeseries
import enesFromTrajFile as eft
import matplotlib.pyplot as plt


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
    """
    """
    g = timeseries.statisticalInefficiency(workArray)
    indices = timeseries.subsampleCorrelatedData(workArray, g)
    workArray_ss = workArray[indices]
    print("%s correlation time is %f" % (label, g))
    return indices, workArray_ss




if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-c","--colvars",
        help="Colvars input file of distance restraint.")
    parser.add_argument("-0", "--traj0",
        help="Data of colvar restraint over trajectory, with restraint OFF.")
    parser.add_argument("-1", "--traj1",
        help="Data of colvar restraint over trajectory, with restraint ON.")

    args = parser.parse_args()
    if None in vars(args).values():
        sys.exit("Error: missing input file(s).")
    print("\nColvars configuration file:\t\t{}".format(args.colvars))
    print("Colvars data file WITHOUT restraint:\t{}".format(args.traj0))
    print("Colvars data file WITH restraint:\t{}\n".format(args.traj1))

    ## fwd work = work required to ADD restraints

    coldict = eft.readColvFile(args.colvars)
    if None in coldict.values():
        sys.exit("Improper colvars configuration file.")

    ### Read data with restraint OFF.
    trajdata_0 = eft.readTrajFile(args.traj0)
    enes_0 = eft.enesFromPositions(trajdata_0, coldict)
    enes_0 = np.asarray(enes_0[1000:]) 

    ### Read data with restraint ON.
    trajdata_1 = eft.readTrajFile(args.traj1)
    enes_1 = eft.enesFromPositions(trajdata_1, coldict)
    enes_1 = -1*np.asarray(enes_1[1001:]) # -1 bc energy of removing restr
    # BAR calls for rev work, and rev work is work to remove restr
    # TODO make skip values (1000) arg parameters

    ### Subsample the work values.
    print("length of original fwd work: {}".format(len(enes_0)))
    print("length of original rev work: {}\n".format(len(enes_1)))
    inds0, work_0_ssd = subsample(enes_0, "no_restrt_energies")
    inds1, work_1_ssd = subsample(enes_1, "restraint_energies")

    ### Apply Bennett acceptance ratio.
    dG, stderr = BAR(work_0_ssd,work_1_ssd) # BAR(fwd, rev)

    ### Output data.
    print("\nSubsampled amounts should be somewhat inversely consistent with correlation times.")
    print("Number of samples (subsampled), fwd work: %d" %len(work_0_ssd))
    print("Number of samples (subsampled), rev work: %d" %len(work_1_ssd))
    print("\nMean work values before subsampling should be somewhat consistent with colvars restraint data.")
    print("mean for fwd work: ",enes_0.mean())
    print("mean for rev work: %f \n" % enes_1.mean())
    print(">>> Free energy, stderr (kcal/mol) to ADD restraints: {0:.4f} +- {1:.4f}\n".format(dG, stderr))
 
    ### Plot (i) restraint data, (ii) work values as line, & subsampled work as dotted line.
    ax1 = plt.subplot(2,2,1) # LHS top
    plt.plot(range(len(enes_0)),enes_0) 
    plt.plot(inds0,work_0_ssd, 'ro', ms=2.0) 

    ax2 = plt.subplot(2,2,3) # LHS bottom
    plt.plot(range(len(enes_1)),enes_1) 
    plt.plot(inds1,work_1_ssd, 'ro', ms=2.0) 

    ax3 = plt.subplot(2,2,2) # RHS top
    plt.plot(range(len(trajdata_0[0])),trajdata_0[0]) 

    ax4 = plt.subplot(2,2,4) # RHS bottom
    plt.plot(range(len(trajdata_1[0])),trajdata_1[0]) 

    ax1.set_title('Adding energy to impose restraint')
    ax2.set_title('Removing energy added from restraint')
    ax3.set_title('Colvars distance w/o restraint')
    ax4.set_title('Colvars distance with restraint')

    ax1.set_xlabel('snapshot')
    ax2.set_xlabel('snapshot')
    ax3.set_xlabel('snapshot')
    ax4.set_xlabel('snapshot')

    ax1.set_ylabel('energy (kcal/mol)')
    ax2.set_ylabel('energy (kcal/mol)')
    ax3.set_ylabel('energy (kcal/mol)')
    ax4.set_ylabel('energy (kcal/mol)')
    plt.tight_layout()
    plt.show()
