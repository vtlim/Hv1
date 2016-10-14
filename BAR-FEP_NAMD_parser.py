#!/usr/bin/env python

### Purpose: Analyze FEP simulations from NAMD.
### Usage: python file.py > file.out

### Assumptions: stepsize is 0.025, 40 windows.
### Note: This script was made in accordance with VMD's ParseFEP plugin:
#    Extensions > Analysis > Analyze FEP simulation
#    with respect to the prob histograms and free energy plots. 
#    Hence the reversed traversal and -1*[] of the reverse FEP data. 
#    If going forward takes x kcal/mol, going backward should take -x kcal/mol.

### Bennett acceptance ratio (BAR) method is used to compute 
#    total free energy difference. See documentation here:
#    https://github.com/choderalab/pymbar/blob/master/pymbar/bar.py

from pymbar import BAR
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob
import re,os


# ------------------------- Variables ---------------------------- #

src = '/work/cluster/limvt/hv1/04_fep/15183_04/F150A/'

# ------------------------- Functions ---------------------------- #


numbers = re.compile(r'(\d+)') # parses a given value
def numericalSort(value):

   """
   Parses some number. 5 would return ['5']. 5.4 would return ['5', '4'].

   """
   parts = numbers.split(value)
   parts[1::2] = map(int, parts[1::2])
   return parts


def cat_fepout(fep_dir, rev=False):

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

   # set name of the summary file
   if rev == True:
      outfile = fep_dir.split('/')[-1]+'_R.fepout'
   else:
      outfile = fep_dir.split('/')[-1]+'_F.fepout'

   # don't write file if already exists
   if os.path.exists(outfile):
      print "File already exists for %s" % (outfile)
      return outfile

   # loop through all *.fepout files and write to the summary file   
   with open(outfile, 'w') as output:
      print outfile
      for fname in fep_file:
         print fname
         with open(fname) as infile:
            output.write(infile.read())

   return outfile


def ParseFEP( fep_file ):

   """
   Parse summary *.fepout files and return relevant data as dictionaries.

   Parameters
   ----------
   fep_file: string. Filename of the summarized results of all 
                     *.fepout results in fep_dir.

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
      if '#STARTING' in l:
         parsing = True

   return dEs_dict, dGs_dict, elecs_dict, vdws_dict, window


def DoBAR(fwds, revs, label):
    """

    BAR to combine fwd and rev data of dGs.
    Here, don't multiply dGs_R by -1 since BAR calls for reverse work value.

    Parameters
    ----------
    fwds: dictionary of forward work values for each window
    revs: dictioary of reverse work values for each window
    label: string label of what it is (only for printing output)

    Returns
    -------
    dgs: 1D list of accumulated list of energy values. Ex. if each step was 2,
       then dgs would be [0,2,4...]
    gsdlist: 1D list of accompanying stdevs to the dgs list
       
    """

    print "\n\nBAR estimator of free energy change:"
    print "\tRel. %s | %s change | uncertainty" % (label, label)
    print "---------------------------------------------------------"
    
    dg_bar = np.zeros([len(fwds)], np.float64)  # allocate storage: dG steps
    gsd_bar = np.zeros([len(fwds)], np.float64) # allocate storage: dG stdev steps
    dgs = np.zeros([len(fwds)], np.float64)     # allocate storage: dG values
    gsdlist = np.zeros([len(fwds)], np.float64) # allocate storage: dG stdev values
    
    for kF, kR in zip(sorted(fwds.keys()), sorted(revs.keys(), reverse=True)):
        dg_bar[kF], gsd_bar[kF] = BAR(fwds[kF],revs[kR])
        print('\t%.4f   %.4f +- %.4f' % (np.sum(dg_bar), dg_bar[kF], gsd_bar[kF]))
    
    # calculate the net dG standard deviation = sqrt[ sum(s_i^2) ]
    gsd = (np.sum(np.power(gsd_bar, 2)))**0.5
    print("Net %s energy difference = %.4f +- %.4f kcal/mol" % (label, np.sum(dg_bar), gsd))
    
    net = 0.
    netsd = 0.
    for i, g in enumerate(dg_bar):
        # accumulate net dGs into running sums (plot this)
        dgs[i] = dg_bar[i] + net
        net = dgs[i]

        # combine the stdevs: s = sqrt(s1^2 + s2^2 + ...)
        gsdlist[i] = ((gsd_bar[i])**2.+(netsd)**2.)**0.5
        netsd = gsdlist[i]


    return dgs, gsdlist


def hist_plot(w_F, w_R, window_F, window_R, title, outfname):

    """
    Plot probability histogram overlap of all windows. Assumes 40 subplots
       from dLambda = 0.025 (1/0.025 = 40 windows).

    Parameters
    ----------

    w_F: dictionary of all the dEs for all windows going forward (key=window number)
    w_R: dictionary of all the dEs for all windows going backward (key=window number)
         Note - based on how FEP calcns are conducted with F and R, need to loop
                over this list in reverse. Aka last of w_R goes with first of w_F.
    window_*: dictionary of start dLambda and stop dLambda per each window. 
            key is the integer lambda window number.
    title: string name of the main title over all 40 windows
    outfname: string name of the image to be saved

    """

    gs = gridspec.GridSpec(8,5)
    idx = 0
    plt.figure()

    for kF, kR in zip(sorted(w_F.keys()), sorted(w_R.keys(), reverse=True)):
        ### set subplot titles based on the dLambda label
#        sbtitle = 'F: %s, R: %s' % (window_F[kF], window_R[kR]) 
        temp = window_F[kF].split(' ') # only uses F. both sets SHOULD match...
        sbtitle = '$\lambda$ = %s to %s' % (temp[1], temp[2])

        plt.subplot(gs[idx])
        plt.hist((-1* w_R[kR] ), bins=50, color='r',histtype='step')
        plt.hist(w_F[kF], bins=50, color='b', histtype='step')
        plt.title(sbtitle, fontsize=6, color='g')
        plt.tick_params(axis='both',labelsize=4)
        idx += 1
    plt.subplots_adjust(top=0.85,hspace=0.85)
    plt.suptitle(title)

    ### super-axis label locations
    plt.text(-20.0, -180.0, '$\Delta$U (kcal/mol)', ha='center') # xlabel
    plt.text(-46.5, 1100.0, 'frequency', va='center', rotation='vertical') # ylabel
    plt.savefig(outfname+'_dE-overlap.eps', format='eps')
    plt.clf()


def dg_plot(dGs_F, dGs_R, window_F, window_R, title, outfname):

    """
    Plot deltaG for all windows. Assumes 40 subplots
       from dLambda = 0.025 (1/0.025 = 40 windows).
    This plot based on equil of 1 ns then FEP data collection 1-5 ns.
       alchEquilSteps 500000, runFEPmin (... nSteps = 2500000 ...)

    Parameters
    ----------
    dGs_F: dictionary of all the dGs for all windows going forward (key=window number)
    dGs_R: dictionary of all the dGs for all windows going backward (key=window number)
         Note - based on how FEP calcns are conducted with F and R, need to loop
                over this list in reverse. Aka last of w_R goes with first of w_F.
    window_*: dictionary of start dLambda and stop dLambda per each window. 
            key is the integer lambda window number.
    title: string name of the main title over all 40 windows
    outfname: string name of the image to be saved

    """

    gs = gridspec.GridSpec(8,5)
    idx = 0
    ns = np.arange(2001)/500.+1
    plt.figure()

    for kF, kR in zip(sorted(dGs_F.keys()), sorted(dGs_R.keys(), reverse=True)):
        ### set subplot titles based on the dLambda label
#        sbtitle = 'F: %s, R: %s' % (window_F[kF], window_R[kR]) 
        temp = window_F[kF].split(' ') # only uses F. both sets SHOULD match...
        sbtitle = '$\lambda$ = %s to %s' % (temp[1], temp[2])

        plt.subplot(gs[idx])
        plt.plot(ns, -1*dGs_R[kR], color='r')
        plt.plot(ns, dGs_F[kF], color='b')
        plt.title(sbtitle, fontsize=6, color='g')
        plt.tick_params(axis='both',labelsize=4)
        idx += 1
    plt.subplots_adjust(top=0.85,hspace=0.85)
    plt.suptitle(title)

    ### super-axis label locations for top=0.85, hspace=0.85
    plt.text(-6.5, -1.2, 'time (ns)', ha='center') # xlabel
    plt.text(-20.0, 9.0, '$\Delta$G (kcal/mol)', va='center', rotation='vertical') # ylabel
    plt.savefig(outfname+'_dG-overlap.eps', format='eps')
    plt.clf()


def gbar_plot(dgs, sds, title, outfname):

    """
    Plot free energy change deltaG over lambda.

    Parameters
    ----------
    dgs: 1D array of dGs
    sds: 1D array of standard deviations corresponding to dgs array.
         If don't have this, just feed function a list of zeroes.
    title: string name of the main title
    outfname: string name of the image to be saved

    """
    ### FOR SOME REASON THE ERROR BARS ARE DISPLAYING HORIZ
    ### EVEN WHEN DEFINING yerr=sds ...

    lambdas = np.linspace(0., 1., len(dgs)) # for x-axis, lambda from 0 to 1
    plt.figure()
    plt.errorbar(lambdas, dgs)
    #plt.errorbar(lambdas, dgs, sds)
    plt.title(title, fontsize=18)
    plt.xlabel("$\lambda$",fontsize=18)
    plt.ylabel("$\Delta$G (kcal/mol)",fontsize=18)
    plt.minorticks_on()
    plt.tick_params(axis='both',width=1.5,length=7,labelsize=16)
    plt.tick_params(which='minor',width=1.0,length=4)
    plt.savefig(outfname+'_summary.eps', format='eps')
    plt.clf()

def pieces_plot(dgs, title, outfname):

    """
    Plot deltaG over lambda as well as the electrostatic and vdW components.

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
    #plt.errorbar(lambdas, dgs, sds)
    plt.title(title, fontsize=18)
    plt.xlabel("$\lambda$",fontsize=18)
    plt.ylabel("energy (kcal/mol)",fontsize=18)
    plt.legend(fancybox=True, loc=2)
    plt.minorticks_on()
    plt.tick_params(axis='both',width=1.5,length=7,labelsize=16)
    plt.tick_params(which='minor',width=1.0,length=4)
    plt.savefig(outfname+'_pieces.eps', format='eps')
    plt.clf()
# ------------------------- Script ---------------------------- #



### Read output data files and summarize results.
fwddir = src+'FEP_F/results'
revdir = src+'FEP_R/results'
fwdfep = cat_fepout(fwddir)
revfep = cat_fepout(revdir, rev=True)
(w_F, dGs_F, elecs_F, vdws_F, window_F) = ParseFEP(fwdfep)
(w_R, dGs_R, elecs_R, vdws_R, window_R) = ParseFEP(revfep)


### Plot probability distributions and energies of fwd and rev.
title = 'FEP Windows: 2GBI in Hv1, F150A (pose 15183_04)'
title += '\nblue = forward (Phe to Ala) | red = reverse (Ala to Phe)'
hist_plot(w_F, w_R, window_F, window_R, title, 'results')
dg_plot(dGs_F, dGs_R, window_F, window_R, title, 'results')

### BAR analysis to combine fwd and rev windows for dG, elec, vdW
alls = np.zeros(shape=(3, len(dGs_F)))
sds = np.zeros(shape=(3, len(dGs_F)))
alls[0], sds[0] = DoBAR(dGs_F, dGs_R, 'dG')
alls[1], sds[1] = DoBAR(elecs_F, elecs_R, 'elec')
alls[2], sds[2] = DoBAR(vdws_F, vdws_R, 'vdW')

## plot BAR summary results
title = 'Free energy change for F150A mutation'
title += '\nof Hv1 (pose 15183_04) with 2GBI bound'
gbar_plot(alls[0], sds[0], title, 'results')
pieces_plot(alls, title, 'z_pieces')
