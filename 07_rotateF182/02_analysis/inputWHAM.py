#!/usr/bin/env python

# Purpose: generate WHAM input file for use in Grossfield WHAM.
#   This is pseudo-command-line script. Change base parameters from US
#   US simulations and for WHAM parameters. Then use command line
#   inputs to specify whether to use subset of data and start / end time (ns). 
# Usage:    python file.py --begin 0 --end 13

# To fix individual traj files (e.g. all +-180 to be around -180):
#  - python, import inputWHAM
#  - inputWHAM.DoChopEndtraj('/path/with/origfile', '/path/with/endfile', 1, 5)

import os
import glob


# ==================================================

### Parameters from US simulations
minZ = -180.0
maxZ = 180.0
spring = 0.1


### WHAM parameters
numbins = 180
tolerance = 0.0001
temp = 300.0
padding = 0

maindir='/pub/limvt/hv1/07_rotateF182/'

os.chdir(maindir)
umbdirs = glob.glob(os.path.join(os.getcwd(),"withF2A","angle*"))
colvarsFile = "colvars.tcl"
trajFile = "npt01.colvars.traj"
chopFile = "npt01.1-5ns.traj"


# ==================================================

def DoChopTraj(trajf, chopf, startns, stopns):
    if startns != 0:
        time1 = (1000*startns/2)+1
    else:
        time1 = startns
    if stopns is not None:
        time2 = (1000*stopns/2)

    with open(trajf,'r') as f:
        lines=f.readlines()

    if os.path.exists(trajf):
    #if os.path.exists(trajf) and not os.path.exists(chopf):
        outf = open(chopf, 'w')
        if stopns is not None: subset = lines[time1:time2]
        else: subset = lines[time1:]
        for i in subset:
            outf.write(str(i))
        outf.close()
        return True
    else: return False


def DoChopEndtraj(trajf, chopf, startns, stopns):
    """
    Not called in the input-wham.py script but wrote to use as a
    standalone function to alter the +- 180 traj files to be all
    around +180 or all around -180 traj files for WHAM.

    Parameters
    ----------
    trajf: string name of input file to be chopped
    chopf: string name of chopped input file
    startns: start ns, assuming 2 fs timestep and colvarsTrajFreq of 1000
    stopns:  stop  ns, assuming 2 fs timestep and colvarsTrajFreq of 1000

    Returns
    -------
    True when successful, False otherwise.

    """
    if startns != 0:
        time1 = (1000*startns/2)+1
    else:
        time1 = startns
    if stopns is not None:
        time2 = (1000*stopns/2)

    with open(trajf,'r') as f:
        lines=f.readlines()

    if os.path.exists(trajf):
    #if os.path.exists(trajf) and not os.path.exists(chopf):
        outf = open(chopf, 'w')
        if stopns is not None: subset = lines[time1:time2]
        else: subset = lines[time1:]

        for line in subset:
            try: value =  float(line.split()[1])
            except ValueError: pass # the step line
            if value < 0: value = value+360  ##### ======== condition to change
            print line.split()[0], value
            outf.write("%s %.14e \n" % (line.split()[0], value))
        outf.close()
        return True
    else: return False


def GenInput(umbdirs, trajFile, chopFile, startns=0, stopns=None):

    ### Output file names
    wham = "wham.in"
    pmf = "pmf.out"
    
    
    ### open and write WHAM input file header
    fname = os.path.join(maindir,'02_analysis','03_wham',wham)
    if not os.path.exists(fname):
        whamf = open(fname, 'w')
    else:
        print("WHAM input file already exists: %s" % fname)
        quit()
    whamf.write("### wham %f %f %d %f %f %f %s %s > wham.out" % (minZ, maxZ, numbins, tolerance, temp, padding, wham, pmf))
    whamf.write("\n### /path/to/timeseries/file loc_win_min spring [correl time] [temp]")
    whamf.write("\n###")
    

    for i, udir in enumerate(umbdirs):
        print udir
        colvf = os.path.join(os.getcwd(),udir,colvarsFile)
        trajf = os.path.join(os.getcwd(),udir,trajFile)
        chopf = os.path.join(os.getcwd(),udir,chopFile)

        if not DoChopTraj(trajf, chopf, startns, stopns):
            print("Error chopping %s." % trajf)

        with open(colvf) as f:
            for line in f:
                if 'centers' in line: cent = line.split()[1]
                if 'forceConstant' in line: spring = line.split()[1]


        whamf.write("\n%s %s %s" % (chopf, cent, spring))
        i += 1
    whamf.close()

# ------------------------- Parse Command Line Inputs ----------------------- #
if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option('-b','--begin',
            help = "Integer start time in nanoseconds. Assumes step*2 / 1e6 = time ns.",
            type = "int",
            dest = 'startns')

    parser.add_option('-e', '--end',
            help = "Integer stop time in nanoseconds. Assumes step*2 / 1e6 = time ns.",
            type = "int",
            dest = 'stopns')

    (opt, args) = parser.parse_args()
    GenInput(umbdirs, trajFile, chopFile, opt.startns, opt.stopns)
