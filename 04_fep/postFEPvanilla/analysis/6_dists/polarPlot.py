
# Purpose: Plot timeseries of a given angle on a polar coordinate plot.
#   Options: (1) scatter plot, (2) histogram of angle.
# Example:
#   python polarPlot.py -i test.dat --hist --time [-f test.fepout] [-j test1.dat test2.dat]

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse

def load_file(filename):
    data = np.loadtxt(filename)
    times = data[:,0]
    angles = data[:,1]
    angles = np.radians(angles)

    return times, angles

def polarTime(filename, fepout=None, skip=0):
    """
    I feel like this is a bit misleading since the radial component is time
      instead of count.

    Reference: https://tinyurl.com/yaj3parp

    """

    # load data
    times, angles = load_file(filename)
    if fepout is not None:
        data_work = np.loadtxt(fepout, comments=['@','&','#'],usecols=6)
        data_work = data_work[0::skip]
        if len(angles) != len(data_work):
            sys.exit("\nERROR: Length of time series and length of fepout data "
                     "do not match.\nCheck that data match, and with a skip "
                     "value of {} for fepout.\n".format(skip))
    # specify colors
    radii = data_work if fepout else times
    cmap   = 'coolwarm' if fepout else 'hsv'

    # generate plot
    ax = plt.subplot(111, projection='polar')
    #c = ax.scatter(angles, times, c=colors, cmap=cmap, alpha=0.55)
    c = ax.scatter(angles, radii, c=times, cmap=cmap, alpha=0.55)
    plt.show()

def polarHist(filename, filelist=[], bins_number=72):
    """
    Reference: https://tinyurl.com/y9wzfo6f
    """

    # load data
    times, angles = load_file(filename)
    # create bins and specify bin widths
    bins = np.linspace(0.0, 2 * np.pi, bins_number + 1)
    width = 2 * np.pi / bins_number
    # histogram data into equally-sized bins and get n=values of bins
    n, _, _ = plt.hist(angles, bins)

    # repeat with filelist if not empty
    nlist = []
    for f in filelist:
        _, angles = load_file(f)
        x, _, _ = plt.hist(angles, bins)
        nlist.append(x)

    # replot the histograms to center bins (e.g., center at 0 instead of start at 0)
    plt.clf()
    ax = plt.subplot(111, projection='polar')
    colors = mpl.cm.Set1(np.linspace(0, 1, 10))
    # IF you have a single histogram & want a diff color, chg index of colors[]
    bars = ax.bar(bins[:bins_number], n, width=width, bottom=0.0,color=colors[0])
    for bar in bars:
        bar.set_alpha(0.5)

    # repeat with filelist if not empty
    for i, f in enumerate(filelist):
        bars = ax.bar(bins[:bins_number], nlist[i], width=width, bottom=0.0, color=colors[i+1])
        for bar in bars:
            bar.set_alpha(0.5)
    # IF you want to change the histogram range, do so here
    ax.set_ylim(0,50)

    for xtick in ax.get_xticklabels():
        xtick.set_fontsize(18)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        help="File with timeseries data. Time in first column,"
                             " measurement in second column. # lines ignored.")
    parser.add_argument("-j", "--jnput", nargs='*',
                        help="MORE files with timeseries data to histogram "
                             "ALONGSIDE the input file specified. Can specify "
                             "multiple files here, separated by space.")
    parser.add_argument("-f", "--fepout",
                        help="NAMD fepout file, from which dE values are used"
                             " to color scatter plot markers.")
    parser.add_argument("-s", "--fepskip", type=int, default=10,
                        help="Take every Nth point of .fepout data. "
                             "Specify 1 to use every data point.")
    parser.add_argument("--hist", action="store_true", default=False,
                        help="")
    parser.add_argument("--time", action="store_true", default=False,
                        help="")

    args = parser.parse_args()
    opt = vars(args)

    # make sure input file is specified AND exists
    if not opt['input'] or not os.path.exists(opt['input']):
        sys.exit("\nERROR: Specify a valid input file.\n")
    # make sure fepout file exists, IF specified
    if opt['fepout'] and not os.path.exists(opt['fepout']):
        sys.exit("\nERROR: Specify valid fepout file.\n")
    # check additional files if specified
    if opt['jnput']:
        for f in opt['jnput']:
            if not os.path.exists(f):
                sys.exit("\nERROR: File '{}' not found.\n".format(f))
    else:
        opt['jnput'] = []

    if args.hist: polarHist(opt['input'], opt['jnput'])
    if args.time: polarTime(opt['input'], opt['fepout'], opt['fepskip'])
