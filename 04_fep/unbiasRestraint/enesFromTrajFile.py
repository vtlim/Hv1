#!/usr/bin/env python

import sys
import csv
import numpy as np

### Written for use of colvars with flat bottom (two half harmonic) potentials
#   defined via lowerBoundary and lowerWallConstant and upper counterparts.
#   To adapt for use by other restraint methods, modify readColvFile dict
#   and enesFromForces function.
### Be careful if reading in energies from file with non-zero width in
#   colvars input file.

def readTrajFile(infile):
    """
    Read in colvar trajectory data from .traj file output by NAMD.

    Parameters
    ----------
    TODO

    Returns
    -------
    TODO

    """

    posList = []
    forList = []

    with open(infile) as fp:
        # consider incorporating delimiters? initial tries didn't work

        rdr = csv.reader((row for row in fp if not row.startswith('#')))
        for row in rdr:
            posList.append(float(row[0].split()[1]))
            try:
                forList.append(float(row[0].split()[2]))
            except IndexError as err:
                pass
    trajData = np.array(list(map(lambda x: np.array(x), [posList, forList])))
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

    enes = list(map(mapfx,trajdata[0]))
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
