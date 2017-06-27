#!/usr/bin/env python

### !!! Python3 Version !!!
### Purpose: Analyze FEP simulations from NAMD using Bennett acceptance ratio.
### Example: python bar4fep.py -d /path/containing/FEPFandFEPR/ -v -p -e 1 -t 10 > output.dat
### Written by: Victoria Lim @ Mobley Lab UCI

### Assumptions: equally distributed fwd and rev windows; 2 fs/step, alchOutFreq = 1000.
### Note: This script was made in accordance with VMD's ParseFEP plugin's
#     dE prob histograms and free energy plots.
#   For rev data, uses reversed traversal and -1*[].
#     If going forward takes x kcal/mol, going backward should take -x kcal/mol.

### Documentation for BAR:
#     https://github.com/choderalab/pymbar/blob/master/pymbar/bar.py

### If ValueError is returned during plot when startStep is specified,
#     make sure (-e and -t) matches range of (startStep to end).
#     E.g., `-e 4 -t 5 -s 1999000`

from pymbar import BAR
from pymbar import timeseries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import re,os,sys,glob
import argparse

# ------------------------- Functions ---------------------------- #


def numericalSort(value):

   """
   Parses some number. 5 would return ['5']. 5.4 would return ['5', '4'].

   """
   numbers = re.compile(r'(\d+)') # parses a given value
   parts = numbers.split(value)
   parts[1::2] = list(map(int, parts[1::2]))
   return parts


def cat_fepout(fep_dir, label, D):

   """
   For a directory containing all lambda windows of forward or
      reverse FEP simulations, combine all results into a single
      file, maintaining numeric order.

   Parameters
   ----------
   fep_dir: string. Full path of direcotry containing FEP results.
   rev: Boolean. True if fep_dir contains reverse FEP results.
                 Default is forward (False).

   Returns
   -------
   outfile: string. Filename of the summarized results of all
                    *.fepout results in fep_dir.

   """
   # get list of all *.fepout file in this fep_dir
   fep_file = sorted(glob.glob(fep_dir+'/*.fepout'), key=numericalSort)
   #outfile = fep_dir.split('/')[-1]+'_{}.fepout'.format(D)
   outfile =  '{}_{}.fepout'.format(label, D)

   # don't write file if already exists
   if os.path.exists(outfile):
      print("!!! WARNING: {} already exists".format(outfile))
      return outfile

   # loop through all *.fepout files and write to the summary file
   with open(outfile, 'w') as output:
      print('Concatenating {}'.format(outfile))
      for fname in fep_file:
         #print(fname)
         with open(fname) as infile:
            output.write(infile.read())

   return outfile


def ParseFEP( fep_file, startStep=None ):

    """
    Parse summary *.fepout files and return relevant data as dictionaries.

    Parameters
    ----------
    fep_file: string. Filename of the summarized results of all
                      *.fepout results in fep_dir.
    startStep: string of integer. Read fepout files starting from this timestep
               (follows the FepEnergy: column)

    Returns
    -------
    dEs_dict: dictionary of all the dE steps for each window (key=window int)
              40 windows = 40 entries. Each key has x values for x num of steps.
    dGs_dict: dictionary of all the dG steps for each window (key=window int)
    window: dictionary of start dLambda and stop dLambda per each window.
            key is the integer lambda window number.

    """
    dEs_dict = {} # dictionary of all the dE steps for each window (key)
    dGs_dict = {} # dictionary of all the dG steps for each window (key)
    elecs_dict = {} # dictionary of all the elec energies for each window (key)
    vdws_dict = {} # dictionary of all the vdw energies for each window (key)
    window = {} # dictionary for the dLambda string labels for each window
    tempDE = [] # temp list to get all dEs in one window
    tempDG = [] # temp list to get all dGs in one window
    tempElec = [] # temp list to get all electrostatic energies of one window
    tempVdw = [] # temp list to get all van der Waals energies of one window
    parsing = False
    i = 0

    # open and get data from fep summary file.
    f = open(fep_file,'r')
    data = f.readlines()
    f.close()

    for line in data:
       l = line.strip().split()

       if '#Free' in l:
          tempDE_arr = np.asarray(tempDE) # convert tempDE to numpy array
          tempDG_arr = np.asarray(tempDG) # convert tempDG to numpy array
          tempElec_arr = np.asarray(tempElec) # convert tempElec to numpy array
          tempVdw_arr = np.asarray(tempVdw) # convert tempVdw to numpy array
          dEs_dict[i] = tempDE_arr # put all the dE values in the dictionary
          dGs_dict[i] = tempDG_arr # put all the dG values in the dictionary
          elecs_dict[i] = tempElec_arr # put all the dE values in the dictionary
          vdws_dict[i] = tempVdw_arr # put all the dG values in the dictionary
          window[i] =  " ".join(l[6:10]) # e.g. grab '[ 0.975 1 ]' join w/space b/t each

          # reset values for the next window of the summary file
          i +=1
          tempDE = []
          tempDG = []
          tempElec = []
          tempVdw = []
          parsing = False

       # append the value in the 'dE' and 'dG' columns of *.fepout file
       if parsing:
          tempDE.append(float(l[6]))
          tempDG.append(float(l[9]))
          tempElec.append(float(l[3])-float(l[2]))
          tempVdw.append(float(l[5])-float(l[4]))

       # turn parsing on at section 'STARTING COLLECTION OF ENSEMBLE AVERAGE'
       if '#STARTING' in l and startStep is None:
          parsing = True
       elif startStep==l[1]: 
          parsing = True

    return dEs_dict, dGs_dict, elecs_dict, vdws_dict, window


def DoBAR(fwds, revs, label, verbose):
    """

    BAR to combine fwd and rev data of dGs.
    Here, don't multiply dGs_R by -1 since BAR calls for reverse work value.

    Parameters
    ----------
    fwds: dictionary of forward work values for each window
    revs: dictionary of reverse work values for each window
    label: string label of what it is (only for printing output)

    Returns
    -------
    dgs: 1D list of accumulated list of energy values. Ex. if each step was 2,
       then dgs would be [0,2,4...]
    gsdlist: 1D list of accompanying stdevs to the dgs list

    """

    fwd_ss = {} # subsampled version of fwds
    rev_ss = {} # subsampled version of revs
    dg_bar = np.zeros([len(fwds)], np.float64)  # allocate storage: dG steps
    gsd_bar = np.zeros([len(fwds)], np.float64) # allocate storage: dG stdev steps
    dgs = np.zeros([len(fwds)], np.float64)     # allocate storage: dG accumulated
    gsdlist = np.zeros([len(fwds)], np.float64) # allocate storage: dG stdev accum


    #corr_time = np.zeros([len(fwds)], np.float64)
    corr_time = {}
    for key, value in fwds.items(): # this notation changes in python3: http://tinyurl.com/j3uq3me
        # compute correlation time
        g = timeseries.statisticalInefficiency(value)
        corr_time[key] = [g]
        # compute indices of UNcorrelated timeseries, then extract those samples
        indices = timeseries.subsampleCorrelatedData(value, g)
        fwd_ss[key] = value[indices]

    for key, value in revs.items(): # this notation changes in python3: http://tinyurl.com/j3uq3me
        # compute correlation time
        g = timeseries.statisticalInefficiency(value)
        corr_time[key].append(g)
        # compute indices of UNcorrelated timeseries, then extract those samples
        indices = timeseries.subsampleCorrelatedData(value, g)
        rev_ss[key] = value[indices]

    bar = {}
    # then apply BAR estimator to get dG for each step
    for kF, kR in zip(sorted(fwd_ss.keys()), sorted(list(rev_ss.keys()), reverse=True)):
        dg_bar[kF], gsd_bar[kF] = BAR(fwd_ss[kF],rev_ss[kR])
        bar[kF] = [ np.sum(dg_bar), dg_bar[kF], gsd_bar[kF] ]

    # calculate the net dG standard deviation = sqrt[ sum(s_i^2) ]
    gsd = (np.sum(np.power(gsd_bar, 2)))**0.5

    net = 0.
    netsd = 0.
    for i, g in enumerate(dg_bar):
        # accumulate net dGs into running sums (plot this)
        dgs[i] = dg_bar[i] + net
        net = dgs[i]
        # combine the stdevs: s = sqrt(s1^2 + s2^2 + ...)
        gsdlist[i] = ((gsd_bar[i])**2.+(netsd)**2.)**0.5
        netsd = gsdlist[i]


    if verbose == True:
        print('\n\n#####---Correlation Times for dG_{}--#####'.format(label))
        print('Window'.rjust(3), 'F'.rjust(5), 'R'.rjust(9))
        for k,v in corr_time.items():
            print("{:3d} {:10.3f} {:10.3f}".format(k, v[0], v[1]) )

        print("\n\n#####---BAR estimator for dG_{}---#####".format(label))
        print('Window'.rjust(3), 'dG'.rjust(5), 'ddG'.rjust(11), "Uncert.".rjust(11))
        print("---------------------------------------------------------")


        for k, v in bar.items():
            str = '{:3d} {:10.4f} {:10.4f} +- {:3.4f}'.format(k, v[0], v[1], v[2])
            print(str)

    print(("\nNet dG_{} energy difference = {:.4f} +- {:.4f} kcal/mol".format(label, np.sum(dg_bar), gsd)))

    return dgs, gsdlist


def hist_plot(w_F, w_R, window_F, window_R, title, outfname):

    """
    Plot probability histogram overlap of all windows. 

    Parameters
    ----------

    w_F: dictionary of all the dEs for all windows going forward (key=window number)
    w_R: dictionary of all the dEs for all windows going backward (key=window number)
         Note - based on how FEP calcns are conducted with F and R, need to loop
                over this list in reverse. Aka last of w_R goes with first of w_F.
    window_*: dictionary of start dLambda and stop dLambda per each window.
            key is the integer lambda window number.
    title: string name of the main title over all windows
    outfname: string name of the image to be saved

    """

    numWins = len(w_F.keys())
    if numWins == 40: gs = gridspec.GridSpec(8,5)
    elif numWins == 20: gs = gridspec.GridSpec(5,4)
    else: print("ERROR: specified number of windows is not currently supported (only 20 or 40)")

    idx = 0
    plt.figure()

    for kF, kR in zip(sorted(w_F.keys()), sorted(list(w_R.keys()), reverse=True)):
        # set subplot titles based on the dLambda label
#        sbtitle = 'F: %s, R: %s' % (window_F[kF], window_R[kR])
        temp = window_F[kF].split(' ') # only uses F. both sets SHOULD match...
        sbtitle = '$\lambda$ = %s to %s' % (temp[1], temp[2])

        # create subplot
        plt.subplot(gs[idx])
        # plot reverse and forward data for this window
        plt.hist((-1* w_R[kR] ), bins=100, color='r',histtype='step')
        plt.hist(w_F[kF], bins=100, color='b', histtype='step')
        # add title and ticks for this window
        plt.title(sbtitle, fontsize=8, color='g')
        plt.tick_params(axis='both',labelsize=6)
        idx += 1
    plt.subplots_adjust(bottom=-1.15,top=1.15,hspace=0.6,wspace=0.3,right=1.3)
    plt.suptitle(title,x=0.7,y=1.3)

    ### super-axis labels. locations are very finicky and likely need adjusting
#    if numWins == 40:
#        plt.text(-20.0, -180.0, '$\Delta$U (kcal/mol)', ha='center') # xlabel
#        plt.text(-46.5, 1100.0, 'frequency', va='center', rotation='vertical') # ylabel
#    if numWins == 20:
#        plt.text(-38.0, -260.0, '$\Delta$U (kcal/mol)', ha='center') # xlabel
#        plt.text(-85.0, 2100.0, 'frequency', va='center', rotation='vertical') # ylabel
    plt.savefig(outfname+'_dE-overlap.eps', format='eps',bbox_inches='tight')
    plt.clf()


def dg_plot(dGs_F, dGs_R, window_F, window_R, eqTime, totTime, title, outfname):

    """
    Plot deltaG for all windows. Based on assumption of 2 fs/step and
        alchOutFreq = 1000 so that each ns has 500 samples.

    Parameters
    ----------
    dGs_F: dictionary of all the dGs for all windows going forward (key=window number)
    dGs_R: dictionary of all the dGs for all windows going backward (key=window number)
         Note - based on how FEP calcns are conducted with F and R, need to loop
                over this list in reverse. Aka last of w_R goes with first of w_F.
    eqTime: float value of how many ns of equilibration per window (based on alchEquilSteps)
    totTime: float value of total sim time per window including eqTime (based on runFEP nSteps)
    window_*: dictionary of start dLambda and stop dLambda per each window.
            key is the integer lambda window number.
    title: string name of the main title over all 40 windows
    outfname: string name of the image to be saved

    """

    numWins = len(dGs_F.keys())
    if numWins == 40: gs = gridspec.GridSpec(8,5)
    elif numWins == 20: gs = gridspec.GridSpec(5,4)
    else: print("ERROR: specified number of windows is not currently supported (only 20 or 40)")

    ### generate ns steps for x-axis
    step=(totTime-eqTime)/((totTime-eqTime)*500.)
    ns=np.arange(eqTime,totTime+step,step)

    idx = 0
    plt.figure()

    for kF, kR in zip(sorted(dGs_F.keys()), sorted(list(dGs_R.keys()), reverse=True)):
        # set subplot titles based on the dLambda label
#        sbtitle = 'F: %s, R: %s' % (window_F[kF], window_R[kR]) # check matching windows
        temp = window_F[kF].split(' ') # only uses F. both sets SHOULD match...
        sbtitle = '$\lambda$ = %s to %s' % (temp[1], temp[2])

        # create subplot
        plt.subplot(gs[idx])
        # plot reverse and forward data for this window
        plt.plot(ns, -1*dGs_R[kR], color='r')
        plt.plot(ns, dGs_F[kF], color='b')
        # add title and ticks for this window
        plt.title(sbtitle, fontsize=8, color='g')
        plt.tick_params(axis='both',labelsize=5)
        # adjust x-axis for this window
        x1,x2,y1,y2 = plt.axis()
        plt.axis((min(ns),max(ns),y1,y2))
        plt.grid()
        idx += 1

    plt.subplots_adjust(bottom=-1.15,top=1.15,hspace=0.6,wspace=0.3,right=1.3)
    plt.suptitle(title,x=0.7,y=1.3)

    ### super-axis labels. locations are very finicky and likely need adjusting
#    if numWins == 40:
#        plt.text(-6.3, -1.2, 'time (ns)', ha='center') # xlabel
#        plt.text(-23.0, 2.0, '$\Delta$G (kcal/mol)', va='center', rotation='vertical') # ylabel
#    if numWins == 20:
#        plt.text(-10.0, -2.5, 'time (ns)', ha='center') # xlabel
#        plt.text(-35.0, 6.5, '$\Delta$G (kcal/mol)', va='center', rotation='vertical') # ylabel
    plt.savefig(outfname+'_dGvTime.eps', format='eps',bbox_inches='tight')
    plt.clf()


def gbar_plot(dgs, sds, title, outfname):

    """
    Plot free energy change deltaG over lambda. 
    Use dG values calculated from BAR.

    Parameters
    ----------
    dgs: 1D array of dGs
    sds: 1D array of standard deviations corresponding to dgs array.
    title: string name of the main title
    outfname: string name of the image to be saved

    """

    lambdas = np.linspace(0., 1., len(dgs)) # for x-axis, lambda from 0 to 1
    plt.figure()
    plt.errorbar(lambdas, dgs, yerr=sds)
    plt.title(title, fontsize=18)
    plt.xlabel("$\lambda$",fontsize=18)
    plt.ylabel("$\Delta$G (kcal/mol)",fontsize=18)
    plt.minorticks_on()
    plt.tick_params(axis='both',width=1.5,length=7,labelsize=16)
    plt.tick_params(which='minor',width=1.0,length=4)
    plt.savefig(outfname+'_summary.eps', format='eps',bbox_inches='tight')
    plt.clf()

def pieces_plot(dgs, title, outfname):

    """
    Plot, as a function of lambda: (1) BAR-calculated deltaG,
       (2) electrostatic component, (3) vdW component.

    Parameters
    ----------
    dgs: array of energies, in order of total dG, elec, vdW
    title: string name of the main title
    outfname: string name of the image to be saved

    """

    lambdas = np.linspace(0., 1., len(dgs[0])) # for x-axis, lambda from 0 to 1
    labels = ['$\Delta$G','electrostatic','van der Waals']
    plt.figure()
    for y, l in zip(dgs, labels):
        plt.plot(lambdas, y, label=l)
    plt.title(title, fontsize=18)
    plt.xlabel("$\lambda$",fontsize=18)
    plt.ylabel("energy (kcal/mol)",fontsize=18)
    plt.legend(fancybox=True, loc=2)
    plt.minorticks_on()
    plt.tick_params(axis='both',width=1.5,length=7,labelsize=16)
    plt.tick_params(which='minor',width=1.0,length=4)
    plt.savefig(outfname+'_decomp.eps', format='eps',bbox_inches='tight')
    plt.clf()
# ------------------------- Script ---------------------------- #

def main(**kwargs):

    def getdata(D, **kwargs):
        ### Read output data files and summarize results.
        src = args.hdir.rstrip('//')
        hdir = src+'/FEP_{}/results'.format(D)
        if os.path.exists(hdir) == True:
            fepout = cat_fepout(hdir, 'results', D)
            (w_D, dGs_D, elecs_D, vdws_D, window_D) = ParseFEP(fepout, args.startStep)
        else:
            raise OSError("No such file or directory '{}'".format(hdir))
        return (w_D, dGs_D, elecs_D, vdws_D, window_D)

    (w_F, dGs_F, elecs_F, vdws_F, window_F) = getdata('F')
    (w_R, dGs_R, elecs_R, vdws_R, window_R) = getdata('R')

    ### BAR analysis to combine fwd and rev windows for dG, elec, vdW
    alls = np.zeros(shape=(3, len(dGs_F))) # actual lists will be shorter bc subsampled
    sds = np.zeros(shape=(3, len(dGs_F)))

    if opt['decomp'] == True:
        alls[2], sds[2] = DoBAR(vdws_F, vdws_R, 'VdW', opt['verbose'])
        alls[1], sds[1] = DoBAR(elecs_F, elecs_R, 'Elec', opt['verbose'])
        alls[0], sds[0] = DoBAR(w_F, w_R, 'Total', opt['verbose'])
    else:
        alls[0], sds[0] = DoBAR(w_F, w_R, 'Total', opt['verbose'])

    ### Plot results.
    if opt['plot'] == True:

        ### Plot probability distributions and energies of fwd and rev.
        print("   Plotting probability distributions...")
        title = 'Energy (dU) Histogram Overlap\nblue = forward | red = reverse'
        hist_plot(w_F, w_R, window_F, window_R, title, opt['outfname'])
        title = 'Free energy (dG) vs. time (ns)\nblue = forward | red = reverse'
        dg_plot(dGs_F, dGs_R, window_F, window_R, float(args.eqTime), float(args.totTime), title, opt['outfname'])

        ### plot BAR summary results
        print("   Plotting free energies...")
        title = "Free energy change over $\lambda$"
        gbar_plot(alls[0], sds[0], title, opt['outfname'])
        if opt['decomp'] == True: pieces_plot(alls, title, opt['outfname'])



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--hdir",
                        help="Location with both FEP_F and FEP_R directories")
    parser.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="Increase output verbosity")
    parser.add_argument("-p", "--plot", action="store_true",
                        help="Generate energy histograms and free energy plots.")
    parser.add_argument("-o", "--outfname", default='plot',
                        help="Base name of saved plots")
    parser.add_argument("-s", "--startStep", default=None,
                        help="Read fepout files starting from this timestep.")
    parser.add_argument("-e", "--eqTime",
                        help="For dG window plot: Nanoseconds of equilibration time per window")
    parser.add_argument("-t", "--totTime",
                        help="For dG window plot: Nanoseconds of total (including equil) time per window")
    parser.add_argument("--decomp", action="store_true", default=False,
                        help="Decompose free energies into electrostatic and vdW components")

    args = parser.parse_args()
    opt = vars(args)
    main(**opt)
