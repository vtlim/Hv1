#!/usr/bin/python

# By: Victoria T. Lim
# If plopping together data from multiple files, Google Sheets can help.
#   Then copy from sheets into vim window.
#   Then replace variable number of spaces with this cmd :%s/ \{2,}/ /g

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
from pymbar import timeseries

# ===========================================


def xyPlot(**kwargs):
    filename = opt['input']
    xlabel = opt['xlabel']
    ylabel = opt['ylabel']
    plttitle = opt['title']
    figname = opt['output']

    with open(filename) as f:
        data = f.read()
    data = data.split('\n')[1:-1] # ignore heading and final newline char *****

    ### Load data for x column.
    leglabel = [(row.split()[0]) for row in data]
    x = [float(row.split()[1]) for row in data]
    numCols = len(data[0].split()[3:])
    numData = len(leglabel)


    ### Load data for y columns.
    y_mat = []
    for i in range(numCols):
       y_mat.append([row.split()[i+3] for row in data]) # ignore label and x columns
#    try:
#        for i in cols[1:]:
#           y_mat.append([row.split()[i] for row in data])
#    except (NameError, IndexError) as err:
#        y_mat.append([row.split()[1] for row in data])
    y_mat = np.array(y_mat)
    y_mat = y_mat.astype(np.float)
    y_vals = y_mat[::2] # get even rows for ddG
    y_stds = y_mat[1::2] # get odd rows for std err

    ### Initialize figure.
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ### Label the figure.
    ax1.set_title(plttitle,fontsize=20)
    ax1.set_xlabel(xlabel,fontsize=18)
    ax1.set_ylabel(ylabel,fontsize=18)
    for xtick in ax1.get_xticklabels():
        xtick.set_fontsize(14)
    for ytick in ax1.get_yticklabels():
        ytick.set_fontsize(14)

    ### Set plot limits.
    axes = plt.gca()
    #axes.set_ylim([-0.5,3])
    axes.set_xlim([min(x)-0.3,max(x)+0.3])


    ### Color the rainbow.
    colors = mpl.cm.gist_rainbow(np.linspace(0, 1, numCols)) # from red to purple
    #colors = mpl.cm.rainbow(np.linspace(0, 0.4, len(y_mat))) # from green to purple
    #colors = mpl.cm.rainbow(np.linspace(0, 0.2, len(y_mat))) # from blue to purple
    #colors = mpl.cm.rainbow(np.linspace(0.4, 1, len(y_mat)))
    markers = ["s","o","8","d","^","x","*","p","v","<","D","+",">","."]

    ### Plot the data.
    # Reference lines
    x.append(min(x)-0.4)
    x.append(max(x)+0.4)
    sortx = np.sort(np.asarray(x))
    ax1.plot(sortx, sortx, color ='gray', linewidth=0.8,label='_nolegend_')
    ax1.plot(sortx, sortx+1, ':', color ='gray', linewidth=0.8, label='_nolegend_')
    ax1.plot(sortx, sortx-1, ':', color ='gray', linewidth=0.8, label='_nolegend_')

    # Points with error bars
    for i, (y, s) in enumerate(zip(y_vals, y_stds)):
        for j, (xi, yi, si) in enumerate(zip(x, y, s)):
            p1 = ax1.errorbar(xi, yi, si, color=colors[j], marker=markers[i],
                 label=leglabel[j], markerfacecolor='None')

    # Legend for mutations without error bars in legend
    patches = [ mpatches.Patch(color=colors[i], label=leglabel[i]) for i in range(numData) ]
    # Custom legend for taut labels
    mark1 = mpl.lines.Line2D([], [], color='black', marker=markers[0],markersize=8)
    mark2 = mpl.lines.Line2D([], [], color='black', marker=markers[1],markersize=8)
    # Generate legend.
    l1 = ax1.legend(handles=patches+[mark1,mark2], labels=leglabel+['taut1','taut2'])
    # Save then show figure.
    plt.savefig(figname, bbox_inches='tight')
    plt.show()


if __name__ == "__main__": 
    import argparse
    parser = argparse.ArgumentParser() 
    parser.add_argument("-i", "--input", 
                        help="Name of the input file.")
    parser.add_argument("-x", "--xlabel",default="",
                        help="Label for x data.")
    parser.add_argument("-y", "--ylabel",default="",
                        help="Label for y data.")
    parser.add_argument("-t", "--title",default="",
                        help="Label for plot title.")
    parser.add_argument("-o", "--output",
                        help="Name of the output figure.")
 
    args = parser.parse_args() 
    opt = vars(args) 
    xyPlot(**opt) 
