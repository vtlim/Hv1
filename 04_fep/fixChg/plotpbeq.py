
"""
Plot electric potential output from CHARMM PBEQ calculation.
These units are labeled as Angstrom for grid and unitCharge/Angstrom for potential.
Verify these with CHARMM configuration file.
By: Victoria Lim

Example: python plotpbeq.py -d /path/to/dir -i -1 -c X

"""

import os
import re
import sys
import pandas as pd
import numpy as np
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt


def get_data(flist):
    """
    Extract data to pandas dataframes from all files in flist.
    All coordinates should be same within all phi files in this directory.

    Parameters
    ----------
    flist - list of phi file names to read from CHARMM PBEQ output

    Returns
    -------
    crds - X Y Z coordinates from 0th 1st 2nd column of first file in flist
    data1 - potentials from 3rd column for all files in flist

    """

    crds = pd.read_csv(flist[0],index_col=None,sep='\s+',header=None,usecols=[0,1,2])
    crds.columns = ["X","Y","Z"]

    templist = []
    for i, f in enumerate(flist):
        if i%10==0: print(i) # status check
        tempdf = pd.read_csv(f,index_col=None,sep='\s+',header=None,usecols=[3])
        templist.append(tempdf)
    data1 = pd.concat(templist, axis=1)
    data1['PHI'] = data1.mean(axis=1)
    df = pd.concat([crds, data1], axis=1)

    return df

def plot_contour(abscissa, ordinate):
    # https://matplotlib.org/examples/pylab_examples/griddata_demo.html
    fig = plt.figure()

    # define grid
    nlevels = 15 # number of contours to draw
    xi = np.linspace(int(min(abscissa)),int(max(abscissa)),1000)
    yi = np.linspace(int(min(ordinate)),int(max(ordinate)),1000)
    zi = mpl.mlab.griddata(abscissa, ordinate, phi, xi, yi, interp='linear')
    cs = plt.contour(xi, yi, zi, nlevels, linewidths=0.5, colors='k')
    cs = plt.contourf(xi, yi, zi, nlevels,
                  vmax=abs(zi).max(), vmin=-abs(zi).max())
    plt.colorbar(label=r"$\Phi$ (e/$\AA$)")
    plt.xlabel('{} ($\AA$)'.format(othr_coord))
    plt.ylabel('Z ($\AA$)')
    plt.title('Frame: {} at {} = 0 $\AA$'.format(frame, coord))
    plt.savefig('{}_{:03d}.png'.format(coord, frame))
#    plt.show()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Plot electrostatic potential from CHARMM PBEQ calculation.')

    parser.add_argument("-d", "--directory", required=True,
        help="Base directory containing .phi files.")
    parser.add_argument("-i", "--frame", required=True, type=int,
        help="Specify frame number for which to plot potential. "
             "To plot average of all frames, use -1")
    parser.add_argument("-p", "--readpickle", default=False,
        help="Look for specified pickle read data from there. "
             "Don't use this flag if reading from .phi files.")
    parser.add_argument("-c","--coordinate", required=True,
        help="Coordinate to view from. Upper case choice of X or Y. "
             "E.g., if specified X, plot will display Y against Z.")

    args = parser.parse_args()
    basedir = args.directory
    frame = args.frame

    # specify coordinate to view from
    if args.coordinate.upper() in set(['X','Y']):
        coord = args.coordinate.upper()
        if coord == 'X': othr_coord = 'Y'
        if coord == 'Y': othr_coord = 'X'
    else:
        sys.exit('Specify a valid coordinate. X Y')

    # generate list of file names
    if frame == -1:
        fs = [x for x in os.listdir(basedir) if x.endswith(".phi")] # find phi files
        fs = ['{}/{}'.format(basedir, f) for f in fs] # add directory to filename
        fs = sorted(fs, key=lambda x: (int(re.sub('\D','',x)),x))
    else:
        fs = "{}/{}.phi".format(basedir,frame)
        if not os.path.isfile(fs):
            sys.exit('ERROR. File not found: {}.'.format(fs))
        fs = [fs]
    print('\n'.join(fs))

    # read in data from files or from pickle
    if not args.readpickle:
        df = get_data(fs)
        pickle.dump(df, open('numpyArrays.pickle','wb'))
    else:
        df = pickle.load(open(args.readpickle, 'rb'))

    # extract data from plane of X=0 or Y=0 (whatever coordinate is specified)
    df_2d = df.loc[df[coord] == 0]
    abscissa = df_2d[othr_coord] # make this generalizable to read in either X or Y
    ordinate = df_2d['Z']
    phi = df_2d['PHI']

    # plot
    plot_contour(abscissa, ordinate)
