#!/usr/bin/python

# By: Victoria T. Lim
# Adapted from my other script: /data12/cmf/limvt/analysis/plotXY.py
# Example:
#   python plotMeasVsWork.py -m meas_R211-CZ_D112-CG.agr -f FEP_R/results/alchemy01.fepout -n 10 -x "every 10k steps" -y "R211:CZ to D112:CG (A)"
# Notes:
#   - Measured data should be plain text file (or similar) with time in x and measurement in y
#   - Can get measured data from calcDist.tcl or from VMD (measure bonds, graph over whole traj, export to Grace or ASCII)

import os, sys
import numpy as np
import random

import matplotlib.pyplot as plt
import matplotlib as mpl

# ===========================================





def plotMeasVsWork(**kwargs):

    in_meas = opt['measurement']
    N = opt['subsample']
    xlabel = opt['xlabel']
    ylabel = opt['ylabel']
    figname = opt['output']

    # if both fepout and work files are specified, exit with error
    if opt['work'] is not None and opt['fepout'] is not None:
        sys.exit("ERROR: Specify work file as text file or fepout file, not both.")
    if opt['work'] is not None:
        in_work = opt['work']
        fromfep = False
    else:
        in_work = opt['fepout']
        fromfep = True

    # open measurement file. ignore lines for comments or grace processing
    data_meas = np.loadtxt(in_meas, comments=['@','&','#'])

    # if work in fepout file, read in 6th column of fepout file (with dE)
    if not fromfep:
        data_work = np.loadtxt(in_work, comments=['@','&','#'])
    else:
        data_work = np.loadtxt(in_work, comments=['@','&','#'],usecols=6)

    # subsample work if specified
    if N is not None:
        data_work = data_work[0::N]
    print(len(data_meas),len(data_work))

    # check that measurement list and work list are the same lengths
    if len(data_meas) != len(data_work):
        sys.exit("ERROR: Number of measurements = {}. Number of work values "
                 "= {}.\nThese should match. Make sure files are consistent "
                 "with each other and a subsample value of {}.".format(
                 len(data_meas),len(data_work), N))

    # check that data do match before plotting
    randint = random.randint(0,len(data_meas))
    yes = {'yes','y', 'ye', ''}
    no = {'no','n'}

    choice = input("At x={}, the measured value is {}, and the work value is "
                   "{}. Are these consistent with each other in the files? "
                   "(Y/N) ".format(randint, data_meas[randint,1],data_work[randint])).lower()
    if choice in yes:
       pass
    elif choice in no:
       sys.exit("Revisit data consistency before plotting.")
    else:
       sys.stdout.write("Please respond with 'yes' or 'no'")


    # convert the x-axis to ns (based on 2 fs step)
#    x = 0.002*np.array(x)


    # plot the data. initialize figure for two subplots.
    fig = plt.figure(figsize=(20,10))
    ax1 = fig.add_subplot(211)
    ax1.plot(data_meas[:,0],data_meas[:,1])
    ax1.set_ylabel(ylabel)

    ax2 = fig.add_subplot(212)
    ax2.plot(data_meas[:,0],data_work)
    ax2.set_ylabel("work (kcal/mol)")
    ax2.set_xlabel(xlabel)

    for ax in (ax1,ax2):
        plt.sca(ax) # set current Axes instance
        plt.grid()
        # find x smallest measurement values
        xcoords = data_work.argsort()[:10]
#        xcoords = data_meas[:,1].argsort()[:10]
        print(" data-sorted x-values with smallest measurements: ",xcoords)
        print("index-sorted x-values with smallest measurements: ",sorted(xcoords))
        for i, xc in enumerate(sorted(xcoords)):
            if abs(xcoords[i-1]-xc) > 1: # don't plot adjacent lines for clarity
                plt.axvline(x=xc,ls='--',c='r',lw=0.8)

    plt.subplots_adjust(hspace=.0)
    plt.savefig('test.eps')
    plt.show()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--measurement", required=True,
                        help="Name of the file with measurements (e.g., "
                        "distances, angles, etc. Two columns: 1st is "
                        "time for x-axis, 2nd is measured value for y-axis.")
    parser.add_argument("-w", "--work",
                        help="Name of file with work values. Should have same "
                        "x-axis as measurements file.")
    parser.add_argument("-f", "--fepout",
                        help="Instead of work file, can specify .fepout file "
                        "from NAMD. That means you should VERIFY that values "
                        "of time/measurements match corresponding work values!")
    parser.add_argument("-n", "--subsample", type=int,
                        help="Take every Nth point of .fepout data.")
    parser.add_argument("-x", "--xlabel",default="",
                        help="Label for x data.")
    parser.add_argument("-y", "--ylabel",default="",
                        help="Label for measurement data.")
    parser.add_argument("-o", "--output",
                        help="Name of the output figure.")


    args = parser.parse_args()
    opt = vars(args)
    plotMeasVsWork(**opt)
