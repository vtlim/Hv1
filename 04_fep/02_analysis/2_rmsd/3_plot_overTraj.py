#!/usr/bin/python

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


# ============ Parameters ===================

pose = '17041_19'
#pose = '17041_13'
#pose = '15183_04'

mut = 'F150A'
way = 'F'
numWins = 40
plotrmsd = True
withLig = True

if plotrmsd:
    filename = "rmsd_endFrames-%s.dat" % way
    figname = "rmsd_endFrames-%s.png" % way
    delimiter = " \t "
    plttitle = "RMSD at end of each %s lambda window,\npose %s, mutation %s" % (way, pose, mut)
    xlabel = "window"
    ylabel = "RMSD ($\AA$)"

    if withLig:
        cols = [1,2] # first and 2nd data columns
        leglabel = ["Hv1 TM backbone", "2GBI"]

    if not withLig:
        numCols = 1 # first n *data* (not time) columns
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



os.chdir('/data12/cmf/limvt/hv1/04_fep/%s/%s/02_analysis/2_rmsd' \
      % (pose, mut))

with open(filename) as f:
    data = f.read()
data = data.split('\n')[1:-1] # -1 gets not a blank line at end

### Generate list for x-axis
x = np.arange(len(data))

### Load data for y columns.
y_mat = []
try:
    for i in cols:
       y_mat.append([row.split(delimiter)[i] for row in data])
except NameError:
    for i in range(1,numCols+1):
       y_mat.append([row.split(delimiter)[i] for row in data])
except IndexError: pass

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
#colors = mpl.cm.rainbow(np.linspace(0, 0.3, n))

### Plot the data.
for color, y in zip(colors, y_mat):
    ax1.scatter(x, y, color=color)
leg = ax1.legend(leglabel,loc=4)

plt.grid()
plt.savefig(figname)
plt.show()


