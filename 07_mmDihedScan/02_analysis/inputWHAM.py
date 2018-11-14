#!/usr/bin/env python

# Purpose: generate WHAM input file for use in Grossfield WHAM.
#   This is pseudo-command-line script. Change base parameters from US
#   US simulations and for WHAM parameters. Then use command line
#   inputs to specify whether to use subset of data and start / end time (ns). 
# Usage:    python file.py --begin 0 --end 13
#   Call this script in the directory with all the subdirectories.
# IMPORTANT: is the spring constant for WHAM the same as for NAMD colvars?
#   if the width is one, yes.
#   if not, take namd forceconstant and divide by width squared.
#   e.g., if width = 0.1 and forceConstant = 0.015, WHAM force constant is 1.5

import os
import glob


# ==================================================

# FOR WHAM INPUT FILE
# http://membrane.urmc.rochester.edu/sites/default/files/wham/doc.pdf

### Parameters from US simulations
minZ = -180.0
maxZ = 180.0
spring = 0.1

### WHAM parameters
numbins = 180
tolerance = 0.0001
temp = 300.0
padding = 0



# ==================================================

def DoChopTraj(trajf, chopf, startns, stopns, translate=False):
    """
    Chops a provided trajectory file based on a given
       start time and end time in nanoseconds. Assuming
       2 fs time step and writing results every 1000 steps.
       Helpful for seeing how PMF evolves over time.

    Parameters
    ----------
    trajf
    translate: Boolean, translate negative angles to positive ones,
        e.g., -90 translates to 270

    Returns
    -------
    True if successful, False otherwise

    """

    if startns != 0:
        time1 = (1000*startns/2)+1
    else:
        time1 = startns
    if stopns is not None:
        time2 = (1000*stopns/2)

    with open(trajf,'r') as f:
        lines=f.readlines()

    # filter out lines to remove commented ones
    filtlines = []
    for line in lines:
        if not line.startswith('#'):
            filtlines.append(line)

    if os.path.exists(trajf):
    #if os.path.exists(trajf) and not os.path.exists(chopf):
        outf = open(chopf, 'w')
        if stopns is not None:
            subset = filtlines[time1:time2]
        else: 
            subset = filtlines[time1:]
        for i in subset:
            value =  float(i.split()[1])
            if translate and value < 0:
                value = value+360  ##### ======== condition to change
            elif not translate and value > 180:
                value = value-360
            outf.write("%s %.14e \n" % (i.split()[0], value))
        outf.close()
        return True
    else:
        print("%s not found in %s" % (trajf, os.getcwd()))
        return False


def GenInput(hdir, trajFile, chopFile, startns=0, stopns=None):
    os.chdir(hdir)

    ### Output file names
    wham = "wham.in"
    pmf = "pmf.out"
    
    
    ### open and write WHAM input file header
    #fname = os.path.join(maindir,'02_analysis','03_wham',wham)
    fname = wham
    if not os.path.exists(os.path.join(hdir,fname)):
        whamf = open(fname, 'w')
    else:
        print("WHAM input file already exists: %s" % fname)
        quit()
    whamf.write("### wham %f %f %d %f %f %f %s %s > wham.out" % (minZ, maxZ, numbins, tolerance, temp, padding, wham, pmf))
    whamf.write("\n### /path/to/timeseries/file loc_win_min spring [correl time] [temp]")
    whamf.write("\n###")
    

    #for i, (dirpath, udir, files) in enumerate(os.walk(hdir)):
    print glob.glob('*/')
    for i, udir in enumerate(glob.glob('*/')):
        print udir
        colvf = os.path.join(os.getcwd(),udir,colvarsFile)
        trajf = os.path.join(os.getcwd(),udir,trajFile)
        chopf = os.path.join(os.getcwd(),udir,chopFile)

        if not DoChopTraj(trajf, chopf, startns, stopns, False):
            print("Error chopping %s." % trajf)

        with open(colvf) as f:
            for line in f:
                if 'centers' in line:
                    cent = float(line.split()[1])
                    if cent > 180: cent = cent-360 #======================================================= change me
                if 'forceConstant' in line: spring = line.split()[1]


        whamf.write("\n%s %f %s" % (chopf, cent, spring))
        i += 1
    whamf.close()

# ------------------------- Parse Command Line Inputs ----------------------- #
if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option('-d','--dir',
            help = "Directory containing all colvar value subdirectories.",
            type = "string",
            dest = 'hdir')

    parser.add_option('-b','--begin',
            help = "Integer start time in nanoseconds. Assumes step*2 / 1e6 = time ns.",
            type = "int",
            dest = 'startns')

    parser.add_option('-e', '--end',
            help = "Integer stop time in nanoseconds. Assumes step*2 / 1e6 = time ns.",
            type = "int",
            dest = 'stopns')

    (opt, args) = parser.parse_args()

    colvarsFile = "colvars.tcl"
    trajFile = "npt01.colvars.traj"
    chopFile = "npt01_1-5ns.traj"

    GenInput(opt.hdir, trajFile, chopFile, opt.startns, opt.stopns)
