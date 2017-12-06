#!/usr/bin/env python

import sys
import csv
import numpy as np

### Written for use of colvars with flat bottom (two half harmonic) potentials
#     defined via lowerBoundary and lowerWallConstant and upper counterparts.
#   To adapt for use by other restraint methods, modify readColvFile dict
#     and enesFromForces function.
### Be careful if reading in energies from file with non-zero width in
#     colvars input file.
### NOT SUITED FOR HANDLING MULTIPLE COLVARS DEFINITIONS, unless all matching

def readTrajFile(infile, col):
    """
    Read in colvar trajectory data from .traj file output by NAMD.
    Index column (0th column) is skipped, and data is read from
    specified column (zero-based indexes).

    Parameters
    ----------
    infile: str, name of input file to read. 0th column is index
      such as time step, all other columns are colvars data.
    col: index of the colvars column to be read.

    Returns
    -------
    TODO

    """

    all_data = np.loadtxt(infile)
    trajData = all_data.T[col] # skip first index column
    return trajData


def readColvFile(infile):
    """
    Read in NAMD colvar input file used to generate .traj data.

    Parameters
    ----------
    infile: string name of NAMD configuration file for colvars module.

    Returns
    -------
    coldict: dictionary of parameters from colvars input file.

    """

    coldict = dict.fromkeys(['width', 'lowerboundary','upperboundary','lowerwallconstant','upperwallconstant'])

    with open(infile) as f:
        inputdata = f.readlines()
    for l in inputdata:
        m = l.split()
        if len(m) == 2 and m[0].lower() in coldict:
            coldict[m[0].lower()] = float(m[1])
    print("Colvars parameters: ",coldict)
    return coldict



def enesFromForces(trajdata, coldict):
    """
    Not developed since outputForce is consistent with position via
    harmonic potential so seems redundant to back calculate
    what I would just put through enesFromPositions anyway.

    """
    if coldict['lowerwallconstant'] != coldict['upperwallconstant']:
        print("This script is not yet adapted to use two difference force constants.")
        return
    pass



def enesFromPositions(trajdata, coldict):

    def harmonic(k,x,x0):
        ene = (1/2.)*k*(x-x0)**2
        return ene

    if coldict['lowerwallconstant'] != coldict['upperwallconstant']:
        print("This script is not yet adapted to use two difference force constants.")
        return
    k = coldict['lowerwallconstant']/coldict['width']/coldict['width']

    mapfx = lambda x: harmonic(k,x,coldict['lowerboundary']) if x < coldict['lowerboundary'] else harmonic(k,x,coldict['upperboundary']) if x > coldict['upperboundary'] else 0

    enes = list(map(mapfx,trajdata))
#    vecfx = np.vectorize(mapfx)  ## returns all zero if the first element fulfills zero condition but works otherwise
#    print(trajdata[0][0:8])
#    print(vecfx(trajdata[0][0:8]))
    return enes



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--infile",
        help="Input the *.traj file with colvars information over the trajectory.")
    parser.add_argument("-c", "--colfile",
        help="Input the input colvars file that was used to generate *.traj file.")
    parser.add_argument("-o", "--outfile",
            help="Name of the output file containing ..... VTL")

    args = parser.parse_args()


    if args.infile is None or args.colfile is None:
        sys.exit("Error: either input file or colvars input file missing.")

    trajdata = readTrajFile(args.infile)
    coldict = readColvFile(args.colfile)
    enes = enesFromPositions(trajdata, coldict)

    import matplotlib.pyplot as plt
    print(enes)
    plt.scatter(range(len(enes)),enes)
    plt.show()
