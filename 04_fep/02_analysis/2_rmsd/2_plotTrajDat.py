#!/usr/bin/python

# Usage: python 2_plotTrajDat.py -f /full/path/with/filename.dat --lig

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import re, argparse
import os

# ===========================================

def main(**kwargs):

    try:
        pose = opt['filename'].split('/')[-5]
    except IndexError:
        print("\nError: Pass in the filename with full path to extract pose information.\n")
        exit()
    mut = opt['filename'].split('/')[-4]

    # how to process the data file
    if opt['chain'] and opt['chain']:
        way = re.search('.byChain\-(.+?).dat',opt['filename']).group(1)
        cols = [1,2,3,4,5,6]
        leglabel = ["Hv1 TM backbone", "S1 backbone","S2 backbone","S3 backbone","S4 backbone","2GBI"]
    elif opt['chain']:
        way = re.search('.byChain\-(.+?).dat',opt['filename']).group(1)
        cols = [1,2,3,4,5] # first five correspond to protein
        leglabel = ["Hv1 TM backbone", "S1 backbone","S2 backbone","S3 backbone","S4 backbone"]
    elif opt['lig']:
        way = re.search('.avgByWin\-(.+?).dat',opt['filename']).group(1)
        cols = [1,2] # first and 2nd data columns
        leglabel = ["Hv1 TM backbone", "2GBI"]
    else:
        way = re.search('.byChain\-(.+?).dat',opt['filename']).group(1)
        #way = re.search('.avgByWin\-(.+?).dat',opt['filename']).group(1)
        numCols = 1 # first n *data* (not time) columns
        leglabel = ["TM backbone"]
    delimiter = " \t "
    figname = os.path.splitext(opt['filename'])[0]+'.png'

    # figure labeling details
    plttitle = "Average RMSD for %s windows\npose %s, mutation %s" % (way, pose, mut)
    xlabel = "window"
    ylabel = "RMSD ($\AA$)"


    with open(opt['filename']) as f:
        data = f.read()
    data = data.split('\n')[1:-1] # -1 gets not a blank line at end

    ### Generate list for x-axis
    x = np.arange(len(data))
    x = x+1 # since no zero window

    ### Load data for y columns.
    y_mat = []
    try:
        for i in cols:
           y_mat.append([row.split(delimiter)[i] for row in data])
    except NameError:
        for i in range(1,numCols+1):
           y_mat.append([row.split(delimiter)[i] for row in data])
    except IndexError: pass

    y_mat = np.array(y_mat)
    ### Initialize figure.
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ### Label the figure.
    ax1.set_title(plttitle,fontsize=20)
    ax1.set_xlabel(xlabel,fontsize=18)
    ax1.set_ylabel(ylabel,fontsize=18)
    for xtick in ax1.get_xticklabels():
        xtick.set_fontsize(16)
    for ytick in ax1.get_yticklabels():
        ytick.set_fontsize(16)

    # set ticks for every Nth window
    plt.xticks(np.arange(min(x)+4, max(x)+5, 5.0))

    ### Set plot limits.
    axes = plt.gca()
    axes.set_ylim([0,6])
    if opt['chain']: axes.set_ylim([0,2])

    ### Color the rainbow.
    n, _ = y_mat.shape
    colors = mpl.cm.tab10(np.linspace(0, 1, 10))
    #colors = mpl.cm.rainbow(np.linspace(0, 0.3, n))

    ### Plot the data.
    for color, y in zip(colors, y_mat):
        ax1.plot(x, y, color=color,lw=1)
    leg = ax1.legend(leglabel,loc=2)


    loc = mpl.ticker.MultipleLocator(base=1)
    ax1.xaxis.set_minor_locator(loc)
    plt.grid(which='both',linewidth=0.3)

    plt.savefig(figname, bbox_inches='tight')
    plt.show()

# ===========================================


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filename",
                        help="Full path filename of *.dat file from RMSD analysis")
    parser.add_argument("--lig", action="store_true", default=False,
                        help="Use flag if ligand RMSD is in *.dat file")
    parser.add_argument("--chain", action="store_true", default=False,
                        help="Use flag if decomposing protein RMSD by chain")

    args = parser.parse_args()
    opt = vars(args)
    main(**opt)
