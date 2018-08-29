
# Example: python joyplotDE.py -e 500000 -f 1000 -l f1_window1.fepout f2_window1.fepout
# Adapted from: https://seaborn.pydata.org/examples/kde_joyplot.html
# By: Victoria Lim

# NOTE: be very careful of histogram option and make sure to manually specify number of bins AND bin range.
#       otherwise relative amounts will be misleading. (TODO: improve this)

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})



def joyplotDE(all_files, inv_files, equil_steps=None, data_freq=None, outname=None, hist=False, nojoy=False):
    """
    all_files : list

    inv_files : list
        List of files whose data should be multiplied by -1 because is reverse work.
    equil_steps : int

    data_freq : int

    outname : string

    hist : Boolean
        If true, generate joyplots by blocky histograms.
        If false, generate joyplots of smoothed KDE plots

    nojoy : Boolean
        If true, plot all curves on same area.
        If false, overlap curves like joyplot format.

    """
    def save_and_return():
        if outname is not None:
            plt.savefig(outname, bbox_inches='tight')

    # Set number of steps to remove for equilibration
    if isinstance(equil_steps, int) and isinstance(data_freq, int):
        lines_remove = int(equil_steps/data_freq) + 2
    else:
        lines_remove = None

    # Load and concatenate files into dataframe
    temp = []
    for i, f in enumerate(all_files):
        print(f)
        d = pd.read_csv(f, delim_whitespace=True, comment='#', header=None, usecols=[6], skiprows=lines_remove)
        if inv_files is not None and f in inv_files:
           d.loc[:,6] *= -1 # the column is named 6 (int) before renaming
        d['label'] = int(i+1)
        temp.append(d)
    df = pd.concat(temp, ignore_index=True)
    df.columns = ['work', 'label']

    # Check out one label's data, without KDE
    #a_values = df.loc[df['label'] == 1, 'work']
    #a_values.plot(kind='hist')
    #plt.show()

    if nojoy:
        colors = ['b','r']
        histoptions = {"histtype": "step", "linewidth": 2, "alpha": 1}
        for i in range(len(all_files)):
            sns.distplot(df.loc[df['label'] == int(i+1)]['work'], color=colors[i], kde=(not hist), hist=hist, hist_kws=histoptions)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel("dE (kcal/mol)",fontsize=14)
        plt.ylabel("normalized probability",fontsize=14)
        save_and_return()
        plt.show()
        return

    # Initialize the FacetGrid object
    pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
    g = sns.FacetGrid(df, row="label", hue="label", aspect=15, size=1., palette=pal)

    if hist:
        # Draw the histograms in the FacetGrid
        g.map(plt.hist, 'work', bins='fd')
    else:
        # Draw the densities in a few steps
        g.map(sns.kdeplot, "work", clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
        g.map(sns.kdeplot, "work", clip_on=False, color="w", lw=2, bw=.2)

    # Add horizontal bottom bar and vertical grid lines
    g.map(plt.axhline, y=0, lw=2, clip_on=False)
    g.map(plt.grid, axis='x', lw=0.5)

    # Define and use a simple function to label the plot in axes coordinates
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color,
                ha="left", va="center", transform=ax.transAxes)
    g.map(label, "work")

    # Set the subplots to overlap
    g.fig.subplots_adjust(hspace=-.25)

    # Remove axes details that don't play well with overlap
    g.set_titles("")
    g.set(yticks=[])
    g.despine(bottom=True, left=True)

    # Show plot
    save_and_return()
#    plt.show()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    # details from simulation
    parser.add_argument("-e", "--equil", type=int, default=500000,
                        help="Number of equilibration steps per window")
    parser.add_argument("-f", "--freq",  type=int, default=1000,
                        help="FEP data was taken every this many steps")
    # details on plots
    parser.add_argument("--hist", action="store_true", default=False,
                        help="Plot histograms instead of KDE distributions")
    parser.add_argument("-o", "--out",
                        help="Name for which to save plot. Leave blank to not save.")
    parser.add_argument("--nojoy", action="store_true", default=False,
                        help="Generate one plot of of all distributions instead"
                             " of separate plots in joyplot format.")
    # input data
    parser.add_argument("-l", "--filelist", nargs='+', # takes 1+ args
                        help="List of individual window files with FEP output data.")
    parser.add_argument("-m", "--multiply", nargs='*', # takes 0+ args
                        help="List of files in filelist that should be multipled"
                             "by negative one due to reverse direction.")

    args = parser.parse_args()
    opt = vars(args)

    joyplotDE(opt['filelist'], opt['multiply'], opt['equil'], opt['freq'], opt['out'], opt['hist'], opt['nojoy'])

