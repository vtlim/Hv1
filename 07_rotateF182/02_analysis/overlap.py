#!/usr/bin/env python

# Purpose: Plot histograms of all umbrella sampling windows altogether; 
#    they should overlap.
# TODO: make command line arguments
# Usage: python file.py

import os
import numpy as np
import matplotlib.pyplot as p

# -------------------- Variables ------------------ #

# end with /
hdir="/pub/limvt/hv1/07_rotateF182/withF2A/"

plotAll = False   # glob and histogram all subdirs

# if plotAll=False, specify range of windows to plot
# full range for f182 is minC=-130, maxC=80
minC=-180
#maxC=-120
maxC=190  # add 'incr' to actual desired stop pt. maxC=150 really means 140.
incr=10

eqt = 0  # equil time to discard. 1000 = 1 ns

trajf = "npt01.colvars.traj"
figname = "fromF2A_all.png"
# ------------------------------------------------- #

os.chdir(hdir)

data = []
for i in range(minC,maxC,incr):
    subdir="angle"+str(i)
    print(subdir)

    tfile=os.path.join(hdir+subdir,trajf)
    if not os.path.isfile(tfile):
        print("%s not found." % (tfile))
        continue

    with open(tfile) as f:
        lines = f.readlines()[eqt:]

        zlist = []
        for line in lines:
            parts = line.split() 
            if not parts[0].startswith("#"): zlist.append(float(parts[1]))
    data.append(zlist)
    
p.figure(figsize=(20,8))
for i, zlist in enumerate(data):
    print(i, len(zlist), min(zlist), max(zlist))
    y,binEdges = np.histogram(zlist,bins=100,range=(-180,180))
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    p.plot(bincenters,y,'-')
p.tick_params(axis='both', which='major', labelsize=18)
p.xlabel("F182 dihedral angle (deg)",fontsize=18)
p.ylabel("count",fontsize=18)

p.savefig(os.path.join(hdir+"../02_analysis/02_overlap/",figname))
p.show()
