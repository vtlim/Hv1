#!/usr/bin/env python

# Purpose: Plot histograms of all umbrella sampling windows altogether; 
#    they should overlap. Reads in colvars output .traj files from NAMD.
# Results setup: single directory with all subdirectories for every value
#    of the collective variable. This script looks for a file named
#    npt01.colvars.traj in each subdirectory.
# Note that python script arguments refer to how to process files. For the
#    min and max values inside the data files, modify np.histogram line in plot.

import os
import numpy as np
import matplotlib.pyplot as p


# ------------------------------------------------- #

def overlap(**kwargs):

    hdir = opt['indir']+'/'
    opt['wmax']=int(opt['wmax'])+int(opt['winc'])
    opt['wmin']=int(opt['wmin'])
    opt['winc']=int(opt['winc'])
    eqt=int(opt['eqt'])
    figname=opt['output']
    
    data = []
    for i in range(opt['wmin'],opt['wmax'],opt['winc']):
    
        tfile=os.path.join(hdir+str(i),'npt01.colvars.traj')
        if not os.path.isfile(tfile):
            print("%s not found." % (tfile))
            continue
        if os.path.getsize(tfile) < 1000:
            print("%s is incomplete." % (tfile))
            continue
    
        with open(tfile) as f:
            lines = f.readlines()[eqt:]
    
            zlist = []
            for line in lines:
                parts = line.split() 
                if not parts[0].startswith("#") and not parts[0].startswith("@"):
                    zlist.append(float(parts[1]))
        data.append(zlist)
        
    p.figure(figsize=(20,8))
    for i, zlist in enumerate(data):
        print(i, len(zlist), min(zlist), max(zlist))
        #y,binEdges = np.histogram(zlist,bins=100,range=(0,360))
        y,binEdges = np.histogram(zlist,bins=100,range=(-180,180))
        bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
        p.plot(bincenters,y,'-')
    p.tick_params(axis='both', which='major', labelsize=18)
    p.xlabel("F182 dihedral angle (deg)",fontsize=18)
    p.ylabel("count",fontsize=18)
    
    p.savefig(os.path.join(hdir,figname))
    p.show()

if __name__ == "__main__": 
    import argparse
    parser = argparse.ArgumentParser() 
    parser.add_argument("-i", "--indir", 
                        help="Name of the input directory with data files.")
    parser.add_argument("-a", "--wmin", 
                        help="Minimum value of all windows. (int)")
    parser.add_argument("-b", "--wmax",
                        help="Maximum value of all windows. (int)")
    parser.add_argument("-c", "--winc",
                        help="Increment between window reference value. (int)")
    parser.add_argument("-e", "--eqt", default=0,
                        help="Equil time to discard. 1000 = 1 ns?") 
    parser.add_argument("-o", "--output", default='overlap.png',
                        help="Name of the output file.")
 
    args = parser.parse_args() 
    opt = vars(args) 
    overlap(**opt) 
