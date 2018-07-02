

# Purpose: Combine lambda windows and plot data for a single direction of a
#             multistage FEP simulation. (either forward or reverse)
# Example: python onedirfep.py -d ../../ -v -p -e 1 -t 5 --revdir -o one_R
# Note:    The code here is a variant of the bar4fep.py script.
#          This script is NOT exact in how it combines windows, but
#             uses a rough approximation by summing the free energy change
#             printed at the bottom of the NAMD fepout files. Main use case
#             of this script is concatenating and plotting.
# Author:  Victoria Lim


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


def cat_fepout(fep_dir, outfile):

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

   # loop through all *.fepout files and write to the summary file
   with open(outfile, 'w') as output:
      print("Number of output files to concatenate for {}: {}".format(outfile, len(fep_file)))
      for fname in fep_file:
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
    elecs_dict: dictionary of the electrostatic component of the dEs
    vdws_dict: dictionary of the vdW component of the dEs
    window: dictionary of start dLambda and stop dLambda per each window.
            key is the integer lambda window number.

    """
    dEs_dict = {} # dictionary of all the dE steps for each window (key)
    dGs_dict = {} # dictionary of all the dG steps for each window (key)
    elecs_dict = {} # dictionary of all the elec energies for each window (key)
    vdws_dict = {} # dictionary of all the vdw energies for each window (key)
    window = {} # dictionary for the dLambda string labels for each window
    namd_net_dg = {} # dictionary for the NAMD dG printed at end of every window
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
          namd_net_dg[i] = l[11]

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

    return dEs_dict, dGs_dict, elecs_dict, vdws_dict, window, namd_net_dg


#    print("\n\n#####---Estimate of single-direction dG_{}---#####".format(label))
#    print('Window'.rjust(3), 'dG'.rjust(5), 'ddG'.rjust(11), "Uncert.".rjust(11))
#    print("---------------------------------------------------------")
#    for k, v in bar.items():
#        str = '{:3d} {:10.4f} {:10.4f} +- {:3.4f}'.format(k, v[0], v[1], v[2])
#        print(str)
#    print(("\nNet dG_{} energy difference = {:.4f} +- {:.4f} kcal/mol".format(label, np.sum(dg_bar), gsd)))


def hist_plot(work_dict, lambda_dict, title, outfname, rev_dir=False):

    """
    Plot probability histogram overlap of all windows.

    Parameters
    ----------

    work_dict: dictionary of all the dEs for all windows (key=window number)
    lambda_dict: dictionary of start dLambda and stop dLambda per each window.
                 key is the integer lambda window number.
    title: string name of the main title over all windows
    outfname: string name of the image to be saved
    rev_dir: Boolean of whether or not the work values are for the reverse
             transformation. Based on how FEP calcns are conducted with fwd
             and rev, the rev data needs to be looped over in reverse. Aka:
             the last data of rev work corresponds with first of fwd work.

    """

    numWins = len(work_dict.keys())
    if numWins == 40: gs = gridspec.GridSpec(8,5)
    elif numWins == 20: gs = gridspec.GridSpec(5,4)
    else: print("ERROR: specified number of windows is not currently supported (only 20 or 40)")

    idx = 0
    plt.figure()

    for k in sorted(list(work_dict.keys()), reverse=rev_dir):
        # set subplot titles based on the dLambda label
        temp = lambda_dict[k].split(' ') # only uses F. both sets SHOULD match...
        sbtitle = '$\lambda$ = %s to %s' % (temp[1], temp[2])

        # create subplot
        plt.subplot(gs[idx])
        # plot data for this window
        if rev_dir:
            plt.hist((-1* work_dict[k] ), bins=100, color='r',histtype='step')
        else:
            plt.hist(work_dict[k], bins=100, color='b', histtype='step')
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


def dg_plot(dg_dict, lambda_dict, eqTime, totTime, title, outfname, rev_dir=False):

    """
    Plot deltaG for all windows. Based on assumption of 2 fs/step and
        alchOutFreq = 1000 so that each ns has 500 samples.

    Parameters
    ----------
    dg_dict: dictionary of all the dGs for all windows (key=window number)
    lambda_dict: dictionary of start dLambda and stop dLambda per each window.
                 key is the integer lambda window number.
    eqTime: float value of how many ns of equilibration per window (based on alchEquilSteps)
    totTime: float value of total sim time per window including eqTime (based on runFEP nSteps)
    title: string name of the main title over all 40 windows
    outfname: string name of the image to be saved
    rev_dir: Boolean of whether or not the work values are for the reverse
             transformation. Based on how FEP calcns are conducted with fwd
             and rev, the rev data needs to be looped over in reverse. Aka:
             the last data of rev work corresponds with first of fwd work.

    """

    numWins = len(dg_dict.keys())
    if numWins == 40: gs = gridspec.GridSpec(8,5)
    elif numWins == 20: gs = gridspec.GridSpec(5,4)
    else: print("ERROR: specified number of windows is not currently supported (only 20 or 40)")

    ### generate ns steps for x-axis
    step=(totTime-eqTime)/((totTime-eqTime)*500.)
    ns=np.arange(eqTime,totTime+step,step)

    idx = 0
    plt.figure()

    for k in sorted(list(dg_dict.keys()), reverse=rev_dir):
        # set subplot titles based on the dLambda label
        temp = lambda_dict[k].split(' ') # only uses F. both sets SHOULD match...
        sbtitle = '$\lambda$ = %s to %s' % (temp[1], temp[2])

        # create subplot
        plt.subplot(gs[idx])
        # plot reverse and forward data for this window
        if rev_dir:
            plt.plot(ns, -1*dg_dict[k], color='r')
        else:
            plt.plot(ns, dg_dict[k], color='b')
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


def gbar_plot(dgs, title, outfname):

    """
    Plot free energy change deltaG over lambda.
    Use dG values printed at end of NAMD output file.

    Parameters
    ----------
    dgs: 1D array of dGs
    title: string name of the main title
    outfname: string name of the image to be saved

    """

    lambdas = np.linspace(0., 1., len(dgs)) # for x-axis, lambda from 0 to 1
    plt.figure()
    plt.plot(lambdas, dgs)
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
    plt.grid()
    plt.savefig(outfname+'_decomp.eps', format='eps',bbox_inches='tight')
    plt.clf()
# ------------------------- Script ---------------------------- #

def main(**kwargs):

    def getdata(D, **kwargs):
        ### Read output data files and summarize results.
        src = args.hdir.rstrip('//')
        hdir = src+'/FEP_{}/results'.format(D)
        outfile =  '{}_{}.fepout'.format('results', D)
        if os.path.exists(outfile):
            print("!!! WARNING: {} already exists".format(outfile))
            (w_D, dGs_D, elecs_D, vdws_D, window_D, namd_D) = ParseFEP(outfile)
        elif os.path.exists(hdir) == True:
            fepout = cat_fepout(hdir, outfile)
            (w_D, dGs_D, elecs_D, vdws_D, window_D, namd_D) = ParseFEP(fepout)
        else:
            raise OSError("No such file or directory '{}'".format(hdir))
        return (w_D, dGs_D, elecs_D, vdws_D, window_D, namd_D)

    if not opt['revdir']:
        (works, dgs, elecs, vdws, windows, namds) = getdata('F')
    else:
        (works, dgs, elecs, vdws, windows, namds) = getdata('R')

    ### Compute list of accumulated energy values (dGs) from the ddGs
    ### Start from the 0 lambdas, correct to fwd dir, sum, fix dir for rev
    ### This section can definitely be improved (todo)
    ddg_list = np.zeros([len(works)], np.float64)     # allocate storage: dG accumulated
    net = 0.
    # i starts at 0 if fwd, i starts at 39 if rev
    for i in sorted(list(namds.keys()), reverse=opt['revdir']):
        # accumulate net dGs into running sums (plot this)
        if opt['revdir']:
            ddg_list[i] = -1.*float(namds[i]) + net
        else:
            ddg_list[i] = float(namds[i]) + net
        net = ddg_list[i]
    if opt['revdir']:
        ddg_list = list(reversed(ddg_list))

    ### Plot results.
    if opt['plot'] == True:

        ### Plot probability distributions and energies of fwd and rev.
        print("   Plotting probability distributions...")
        title = 'Energy (dU) Histogram Overlap'
        hist_plot(works, windows, title, opt['outfname'], opt['revdir'])
        title = 'Free energy (dG) vs. time (ns)'
        dg_plot(dgs, windows, float(args.eqTime), float(args.totTime), title, opt['outfname'], opt['revdir'])

        ### plot BAR summary results
        print("   Plotting free energies...")
        title = "Free energy change over $\lambda$"
        gbar_plot(ddg_list, title, opt['outfname'])
#        if opt['decomp'] == True: pieces_plot(alls, title, opt['outfname'])



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--hdir",
                        help="Location with FEP directories")
    parser.add_argument("-r", "--revdir", action="store_true", default=False,
                        help="Data is from reverse transformation.")
    parser.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="Increase output verbosity")
    parser.add_argument("-p", "--plot", action="store_true",
                        help="Generate energy histograms and free energy plots.")
    parser.add_argument("-o", "--outfname", default='oneplot',
                        help="Base name of saved plots")
    parser.add_argument("-e", "--eqTime",
                        help="For dG window plot: Nanoseconds of equilibration time per window")
    parser.add_argument("-t", "--totTime",
                        help="For dG window plot: Nanoseconds of total (including equil) time per window")
#    parser.add_argument("--decomp", action="store_true", default=False,
#                        help="Decompose free energies into electrostatic and vdW components")

    args = parser.parse_args()
    opt = vars(args)
    main(**opt)
