#!/usr/bin/python

# By: Victoria T. Lim
# If plopping together data from multiple files, Google Sheets can help.
#   Then copy from sheets into vim window.
#   Then replace variable number of spaces with this cmd :%s/ \{2,}/ /g

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pymbar import timeseries

# ===========================================



def subSample(x, y_mat):
    """
    Parameters
    ----------
    Returns
    -------
    """
    #x_mat=np.empty([len(y_mat),len(x)])
    x_mat = []
    z_mat = [] # the subsampled y_mat
    for i, y in enumerate(y_mat):
        # Compute correlation times.
        g = timeseries.statisticalInefficiency(y)
        indices = timeseries.subsampleCorrelatedData(y, g)
        # Subsample data.
        y_sub = y[indices]
        x_sub = x[indices]
        z_mat.append(y_sub)
        x_mat.append(x_sub)
        print("\nLength of original timeseries data: %d\nLength of subsampled\
 timeseries data: %d\n" % (len(y), len(y_sub)) )
    return x_mat, z_mat


def runningMean(y_mat, N):
    # http://stackoverflow.com/questions/13728392/moving-average-or-running-mean
    z_mat = [] # the subsampled y_mat
    for i, y in enumerate(y_mat):
        y_mean = np.convolve(y, np.ones((N,))/N,mode='valid')
        z_mat.append(y_mean)
    return z_mat




#    if legend: leg = ax1.legend(leglabel,loc='center right',bbox_to_anchor=(1.4,0.5))  # legend on outside right
#    plt.savefig(figname, bbox_extra_artists=(leg,), bbox_inches='tight') # bbox is xy where origin is left bottom


def xyPlot(**kwargs):

    def formatFig(ax1, plt):

        if publish:
            tSize = 12
            xSize = 8
            ySize = 8
            kSize = 8
            fig.set_size_inches(3.37,1.7)
            if legend:
#                leg = ax1.legend(leglabel,loc='upper center',fontsize=8,bbox_to_anchor=(0.50,1.65),ncol=2) # all PMFs together
                leg = ax1.legend(leglabel,loc='upper center',fontsize=8,bbox_to_anchor=(0.50,1.30),ncol=2)
        else:
            tSize = 20
            xSize = 18
            ySize = 18
            kSize = 16
            if legend:
                leg = ax1.legend(leglabel,loc='upper left')

        ### Label the figure.
        ax1.set_title(plttitle,fontsize=tSize)
        ax1.set_xlabel(xlabel,fontsize=xSize)
        ax1.set_ylabel(ylabel,fontsize=ySize)
        for xtick in ax1.get_xticklabels():
            xtick.set_fontsize(kSize)
        for ytick in ax1.get_yticklabels():
            ytick.set_fontsize(kSize)

        ### Save figure.
        if publish:
            plt.savefig(figname, bbox_inches='tight',dpi=300)
        else:
            plt.savefig(figname, bbox_inches='tight')

        plt.show()

    def readFile(filename):
        """
        """
        with open(filename) as f:
            data = f.read()
        data = data.split('\n')[1:-1] # accounts for heading and final newline char *****
        return data

    def parseColumns(data):
        y_mat = []
        try:
            for i in cols:
               y_mat.append([row.split(delimiter)[i] for row in data])
        except (NameError, IndexError) as err:
            y_mat.append([row.split(delimiter)[1] for row in data])
        y_mat = np.array(y_mat)
        y_mat = y_mat.astype(np.float)
        return y_mat


    filename = opt['input']
    delimiter = opt['delimiter']
    uncertf = opt['uncert']
    doSubsample =  opt['subsample']
    xlabel = opt['xlabel']
    ylabel = opt['ylabel']
    plttitle = opt['title']
    figname = opt['output']
    publish =  opt['publish']
    if opt['mean'] != 0:
        runMean = True
        runLength = int(opt['mean'])
    else:
        runMean = False
    if opt['columns'] is not None:
        cols = list(map(int,opt['columns'].split(',')))
    if opt['legend'] is not None:
        leglabel = opt['legend'].split(';')
        legend=True

    ### Read in data.
    data = readFile(filename)
    if uncertf is not None:
        uncerts = readFile(uncertf)

    ### Load data for x column.
    x = [float(row.split(delimiter)[0]) for row in data]

#    # Convert the x-axis to ns (based on 2 fs step)
#    x = 0.002*np.array(x)

    ### Load data for y columns.
    y_mat = parseColumns(data)
    if uncertf is not None:
        u_mat = parseColumns(uncerts)

    ### subsample data (may not want to if not timeseries data!)
    if doSubsample: x_mat, y_mat = subSample(x,y_mat)
    if runMean:
        y_mat = runningMean(y_mat,runLength)
        x = 0.002*np.asarray(range(len(y_mat[0])),dtype=np.float32) # now x is approximate, not exactly matching with y

    ### Initialize figure.
    fig = plt.figure()
    ax1 = fig.add_subplot(111)


    ### Set plot limits.
    axes = plt.gca()
    axes.set_ylim([-0.1,3])
#    axes.set_xlim([min(x)-10,max(x)+10])
    axes.set_xlim([min(x)-0.2,max(x)+0.2])


    ### Color the rainbow.
#    colors = mpl.cm.Set1(np.linspace(0, 1, 8)) # qualitative
    colors = mpl.cm.tab20(np.linspace(0, 1, len(y_mat)+5)) # qualitative
    #colors = mpl.cm.rainbow(np.linspace(0, 0.4, len(y_mat))) # from green to purple
    #colors = mpl.cm.rainbow(np.linspace(0, 0.2, len(y_mat))) # from blue to purple
    #colors = mpl.cm.rainbow(np.linspace(0.4, 1, len(y_mat)))

    # manually editing colors
#    colors=colors[:2]
#    colors=[colors[1],colors[3]]

    ### Plot the data.
    if doSubsample:
        for color, x, y in zip(colors, x_mat, y_mat):
            ax1.plot(x, y, color=color)
    elif uncertf is not None:
        for color, y, u in zip(colors, y_mat, u_mat):
#            ax1.errorbar(x, y, yerr=u, color=color)
            ax1.errorbar(x, y, yerr=u, capsize=0.8,lw=0.8,color=color)
    else:
        for color, y in zip(colors, y_mat):
            print(len(x),len(y))
            ax1.plot(x, y, color=color)
#            ax1.plot(x, y, lw=0.8, color=color) # thinner line

    ### Custom text on plot
#    ax1.text(2,11,"A",fontsize=10)

    formatFig(ax1,plt)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        help="Name of the input file. First line is assumed "
                        + "to be some heading line and is NOT read in.")
    parser.add_argument("-d", "--delimiter",
                        help="Put in quotes the delimiter separating columns.")
    parser.add_argument("-c", "--columns",default=None,
                        help="Specify particular data columns to plot.\
                        Separate values with commas.\
                        0th column is x, so don't specify 0.\
                        If not specified, will only plot first data column.")
    parser.add_argument("-u", "--uncert",default=None,
                        help="Name of the file with corresponding uncertainties"
                        + ". Not compatible with running means or subsampling.")
    parser.add_argument("-m", "--mean", default=0,
                        help="If not default=0, take running means over the "
                        + "specified number of data points for each column.")
    parser.add_argument("-s", "--subsample", action="store_true",default=False,
                        help="Subsample y data based on correlation times.")
    parser.add_argument("-x", "--xlabel",default="",
                        help="Label for x data.")
    parser.add_argument("-y", "--ylabel",default="",
                        help="Label for y data.")
    parser.add_argument("-t", "--title",default="",
                        help="Label for plot title.")
    parser.add_argument("-l", "--legend",default=None,
                        help="Add legend to data. Format input as data1;data2;...")
    parser.add_argument("-o", "--output",
                        help="Name of the output figure.")
    parser.add_argument("--publish", action="store_true",default=False,
                        help="Reduce figure/font sizes for publications.")

    args = parser.parse_args()
    opt = vars(args)
    xyPlot(**opt)
