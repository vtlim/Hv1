#!/usr/bin/python

### Purpose: process energy summary file
### Usage: python file.py > output.dat
#python dihedParser.py -d /path/containing/QMandMM/ -q output.dat -m minimize.log -k 2000 --show --save


import os, glob, re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse


# --------------------------------------------- #

def numericalSort(value):

   """
   Parses some number. 5 would return ['5']. 5.4 would return ['5', '4'].
   Only used for sorting in cat_* functions.

   """
   numbers = re.compile(r'(\d+)') # parses a given value
   parts = numbers.split(value)
   parts[1::2] = list(map(int, parts[1::2]))

   return parts



def ReadSumFile(filename):
    angs = []
    enes = []
    with open(filename) as ff:
        for line in ff:
            parts = line.split()
            angs.append(float(parts[0]))
            enes.append(float(parts[1]))
    angs = np.asarray(angs)
    enes = np.asarray(enes)
    return angs, enes


def cat_QM(ddir, fname, qOther='dihed-qm', theory='mp2-631Gd'):

    """
    For a directory containing subdirectories of all angles,
       get the final energy of the angle from the output fname.
 
    Parameters
    ----------
    ddir: string | full path of directory containing *both* QM and MM angle directories.
    fname: string | name of the output file. assumed same for all angle jobs.
    qOther: string | name of the QM dir within ddir. In case is not regularly named.
    theory: string | part of the path with the level of theory identifier. 
 
    Returns
    -------
    angs: numpy array of reference dihedral angle of scan
    enes: numpy array with final energy of QM optimization (if completed)
 
    """

    qmdir = ddir+'/'+qOther

    # Get list of all angles' output files
    # One * for angle and one * for level of theory
    qmfiles = sorted(glob.glob(qmdir+'/*/'+theory+'/'+fname), key=numericalSort)
    outfile = 'summary-qm.dat'
    angs = np.zeros(len(qmfiles))
    enes = np.zeros(len(qmfiles))

    # don't write file if already exists
    if os.path.exists(outfile):
       print("!!! WARNING: {} already exists. Reading from summary file.".format(outfile))
       angs, enes = ReadSumFile(outfile)
       return angs, enes

    # open the output file for summarized results
    with open(outfile, 'w') as output:
        # loop over all angle output files
        for i, filename in enumerate(qmfiles):
            with open(filename) as fname:
                # for this file, look for the line with the final energy
                for line in fname:
                    if "Final energy" in line:
                        angle = filename.split(qOther)[1].split('/')[1]
                        energy = float(re.sub(' +', ' ', line).split(' ')[3])

                        angs[i] = angle
                        enes[i] = energy
                        output.write('{}\t{}\n'.format(angle, energy))
    return angs, enes




def cat_MM(ddir, fname):

    """
    For a directory containing subdirectories of all angles,
       get the final energy of the angle from the output fname.
 
    Parameters
    ----------
    ddir: string | full path of directory containing *both* QM and MM angle directories.
    fname: string | name of the output file. assumed same for all angle jobs.
 
    Returns
    -------
    angs: numpy array of reference (not actual) dihedral angle of scan
    enes: numpy array with potential energy term from NAMD log files
 
    """
    # Get list of all angles' output files
    mmdir = ddir+'/dihed-mm'
    mmfiles = sorted(glob.glob(mmdir+'/*/'+fname), key=numericalSort)
    outfile = 'summary-mm.dat'
    angs = np.zeros(len(mmfiles))
    enes = np.zeros(len(mmfiles))

    # don't write file if already exists
    if os.path.exists(outfile):
       print("!!! WARNING: {} already exists. Reading from summary file.".format(outfile))
       angs, enes = ReadSumFile(outfile)
       return angs, enes

    # open the output file for summarized results
    with open(outfile, 'w') as output:
        # loop over all angle output files
        for i, filename in enumerate(mmfiles):
            # for this file, look for the line with the final energy
            for line in reversed(open(filename).readlines()):
                if line.startswith('ENERGY:'):
                    angle = filename.split('dihed-mm')[1].split('/')[1]
                    energy = float(re.sub(' +', ' ', line).split(' ')[13])

                    angs[i] = angle
                    enes[i] = energy
                    output.write('{}\t{}\n'.format(angle, energy))
                    break
    return angs, enes



def SubtractRestraintE(filename, fConst):
    """

    Subtract restraint energy for MM energies. Using harmonic (not cos)
    potential based on setup in the extra bonds file. 
                 k * (x_actual - x_reference) ^ 2

    Parameters
    ----------
    filename: name of the file with actual, measured dihedral angles from minimization
    fConst: value of the spring constant k declared in NAMD extra bonds file

    Returns
    -------
    xact
    rpe: numpy array with restraint potential energies for all dihedral angles in filename

    """
    xref = []   # reference dihedral angles
    xact = []   # actual dihedral angles

    with open(filename) as ff:
        for line in ff:
            parts = line.split()
#            i = int(parts[0])/5  # get index for list. assumes ALL angles are present & in order
            xref.append( parts[0] )
            xact.append( parts[1] )

    # Convert to numpy arrays
    xref = np.asarray(xref, dtype=np.int32)
    xact = np.asarray(xact, dtype=np.float32)
    
    # If the angle is negative, add 360.
    for i,ang in enumerate(xact):
        if ang < 0: xact[i] = ang+360
    
    # for the first angle, make ~0 if it's ~360
    if abs(xact[0]) > 5:
        xact[0] = xact[0] - 360
    
    # Get the restraint energies.
    rpe = float(fConst)*np.subtract(xact,xref)**2.
    return xact, rpe


def plotDihedScan(x, y, pdict, toSave, toShow, x2=None, y2=None):
    """
    Generate scatter plot given lists of xy data. Optionally, can 
      plot two sets of data on the same plot.

    Parameters
    ----------
    x: list of floats
    y: list of floats
    pdict: dictionary with plot details. should contain keys for
      'title','figname'. If plotting 2, should also contain
      'color1','color2','label1','label2'.
    toSave: Boolean, save figure
    toShow: Boolean, show figure
    x2: second list of x floats, optional 
    y2: second list of y floats, optional
    
    """

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    xlabel = "angle (degrees)"
    ylabel = "energy (kcal/mol)"

    ### Label the figure, larger font
    ax1.set_title(pdict['title'],fontsize=20)
    ax1.set_xlabel(xlabel,fontsize=18)
    ax1.set_ylabel(ylabel,fontsize=18)

    ### Increase font size of tick labels
    #ax1.set_xticklabels(xticks,fontsize=14)
    for ytick in ax1.get_yticklabels():
        ytick.set_fontsize(14)
    for xtick in ax1.get_xticklabels():
        xtick.set_fontsize(14)

    if x2 is None and y2 is None:
        ax1.scatter(x, y,edgecolors='none')

    elif x2 is not None and y2 is not None:
        ax1.scatter(x, y, edgecolors='none',  color=pdict['color1'],label=pdict['label1'])
        ax1.scatter(x2, y2, edgecolors='none',color=pdict['color2'],label=pdict['label2'])

    else:
        print("Either none or both x2 and y2 should be defined.")
        return

    plt.grid()
    plt.grid(which='both', color='0.65',linestyle='-')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    if toSave == True: plt.savefig(pdict['figname'],bbox_inches='tight')
    if toShow == True: plt.show()
    plt.ioff()
    plt.close()
    

# ------------------------- Script ---------------------------- #

def main(**kwargs):

    ###  QM
    ang_qm, ene_qm = cat_QM(opt['ddir'], opt['qfile'], 'dihed-180',opt['theory'])
    #ang_qm, ene_qm = cat_QM(opt['ddir'], opt['qfile'],opt['theory'])

    indices = np.nonzero(ene_qm)
    ang_qm = ang_qm[indices]
    ene_qm = ene_qm[indices]

    minE = min(ene_qm)
    rel_ene_qm = [627.5095*(i - minE) for i in ene_qm] # convert Hartrees -> kcal/mol

    pdict = {}
    pdict['title'] = "Dihedral Scan for 2GBI Tautomer #2 - QM"
    pdict['figname'] = "plot_relDihed-qm.png"
    plotDihedScan( ang_qm, rel_ene_qm, pdict, opt['save'], opt['show'] )


    ###  MM
    if not opt['qOnly']:
        ang_mm, ene_mm = cat_MM(opt['ddir'], opt['mfile'])
        xact, rpe = SubtractRestraintE(opt['ddir']+'/dihed-mm/diheds-from-coor.dat',opt['fConst'])
    
        ### Subtract actual energies minus restraint energies.
        rs_ene_mm = np.subtract(ene_mm,rpe)
    
        ### Take relative energies from minimum.
        minE = min(ene_mm)
        rel_ene_mm = [i - minE for i in ene_mm]
        minE = min(rs_ene_mm)
        rel0_ene_mm = [i - minE for i in rs_ene_mm]

        ### Plot MM results.
        pdict['title'] = "Dihedral Scan for 2GBI Tautomer #2 - MM"
        pdict['figname'] = "plot_relDihed-mm.png"
        plotDihedScan( ang_mm, rel_ene_mm, pdict, opt['save'], opt['show'] )

        ### Plot both QM and MM results.
        pdict['title'] = "Dihedral Scan for 2GBI Tautomer #2"
        pdict['figname'] = "plot_relDihed.png"
        pdict['label1'] = "MM (NAMD, CGenFF)"      # MM line label
        #pdict['label2'] = "QM (Psi4, MP2/6-31G*)"  # QM line label
        pdict['label2'] = "QM (Psi4, MP2/def2-tzvp)"  # QM line label
        pdict['color1'] = 'b'      # MM line color
        pdict['color2'] = 'r'      # QM line color
        plotDihedScan( ang_mm, rel_ene_mm, pdict, opt['save'], opt['show'], ang_qm, rel_ene_qm)

    ### Write out results.
    if not os.path.exists("summary.dat"):
        with open("summary.dat",'w') as writeout:

            writeout.write("\n\tRESULTS FOR QM DIHEDRAL SCAN")
            writeout.write("\n\t----------------------------")
            writeout.write("\nangle_ref\ttotE(Har)\t\trelE(kc/mol)")
            for i in range(len(rel_ene_qm)):
                writeout.write("\n\t%.1f\t\t%.3f\t%.3f" % (ang_qm[i], ene_qm[i], rel_ene_qm[i]))

            if not opt['qOnly']:
                writeout.write("\n\n\tRESULTS FOR MM DIHEDRAL SCAN")
                writeout.write("\n\t----------------------------")
                writeout.write("\nangle_ref\tangle_act\ttotalE\t\trestrE\t(tot-restr)\torigTotE (rel)")
                for i in range(len(rel_ene_mm)):
                    writeout.write("\n\t%.2f\t%f\t%.3f\t%.3f\t%.3f\t%.3f" % (ang_mm[i], xact[i], ene_mm[i], rpe[i], rs_ene_mm[i],rel_ene_mm[i]))
    else:
        print("!!! WARNING: {} already exists. Skip writing summary results.".format('summary.dat'))




#    ### Plot total, restraint, and tot-restr energies for MM scan.
#    fig = plt.figure()
#    ax1 = fig.add_subplot(111)
#    plttitle = "Energies of MM dihedral Scan"
#    FormatPlot(plttitle)
#    figname = "plot_totRestrEnes.png"
#    ax1.plot(x, rpe, label='harmonic restraint E')
#    ax1.plot(x, enes, label='total original E')
#    ax1.plot(x, fins, label='total - restraint E')
#    plt.legend(loc='center right')
#    if opt['save'] == True: plt.savefig(figname,bbox_inches='tight')
#    if opt['show'] == True: plt.show()
#    plt.ioff()
#    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--ddir",
                        help="Location with both QM and MM directories")

    # QM-related arguments
    parser.add_argument("-q", "--qfile",
                        help="Name of the QM output file (assuming same name for all angles)")
    parser.add_argument("-t", "--theory", default='mp2-631Gd',
                        help="Part of the path with the level of theory identifier.")
    parser.add_argument("--qOnly", action="store_true", default=False,
                        help="Only process QM results.")

    # MM-related arguments
    parser.add_argument("-m", "--mfile",
                        help="Name of the MM output file (assuming same name for all angles")
    parser.add_argument("-k", "--fConst",
                        help="Numeric value of the force constant used to restrain dihedral in MM.\
             Is NOT currently used, since dihedral scan looks more reasonable without that value.")

    # Plot-related arguments
    parser.add_argument("--show", action="store_true", default=False,
                        help="Display all plots generated.")
    parser.add_argument("--save", action="store_true", default=False,
                        help="Save all plots.")


    args = parser.parse_args()
    opt = vars(args)
    main(**opt)

