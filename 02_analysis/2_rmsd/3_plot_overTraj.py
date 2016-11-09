#!/usr/bin/python

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


# ============ Parameters ===================

#pose = '17041_19'
#pose = '17041_13'
pose = '15183_04'

mut = 'F150A-noGBI'
numWins = 20
plotrmsd = True

if plotrmsd:
    filename = "rmsd_endFrames.dat"
    figname = "rmsd_endFrames-%s.png" % pose
    figname = "rmsd_endFrames.png"
    delimiter = " \t "
    numCols = 1 # first n *data* (not time) columns
#    cols = [1,3] # first and 3 data columns
    xlabel = "window"
    ylabel = "RMSD ($\AA$)"
    plttitle = "RMSD of Hv1 Backbone (no 2GBI), pose %s" % (pose)
    leglabel = ["TM backbone"]

if not plotrmsd:
    filename = "hv1+gbi_contacts.dat"
    delimiter = "  "
    numCols = 5
    xlabel = "time (ns)"
    ylabel = "Distance ($\AA$)"
    plttitle = "Hv1 Contacts with 2GBI, pose %s" % pose
    leglabel = ["F150-benzo","R211-guan","D112-imid","S181-imid","R211-imid"]
    figname = "plot_contacts_%s.png" % pose

# ===========================================



os.chdir('/data12/cmf/limvt/hv1/04_fep/%s/%dwindows/%s/02_analysis/2_rmsd' \
      % (pose, numWins, mut))

with open(filename) as f:
    data = f.read()
data = data.split('\n')[1:-1] # -1 gets not a blank line at end

### Generate list for x-axis
x = np.arange(len(data))+1

### Load data for y columns.
y_mat = []
try:
    for i in cols:
       y_mat.append([row.split(delimiter)[i] for row in data])
except NameError:
    for i in range(1,numCols+1):
       y_mat.append([row.split(delimiter)[i] for row in data])

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

# set ticks for every other window
plt.xticks(np.arange(min(x), max(x)+1, 2.0))

### Set plot limits.
#axes = plt.gca()
#axes.set_ylim([0,6])

### Color the rainbow.
n, _ = y_mat.shape
colors = mpl.cm.rainbow(np.linspace(0, 1, n))
#colors = mpl.cm.rainbow(np.linspace(0, 0.4, n))

### Plot the data.
for color, y in zip(colors, y_mat):
    ax1.scatter(x, y, color=color)
leg = ax1.legend(leglabel,loc=4)

plt.savefig(figname)
plt.show()


