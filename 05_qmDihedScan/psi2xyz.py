#!/usr/bin/python

### Purpose: Grab last set of optimized coordinates from Psi4 output file and 
###    write to XYZ file format.
### Usage: python file.py -i psi4output.dat -o label.xyz -t pdb

import argparse

def main(**kwargs):

    keeper = False
    keeplines = []
    for line in reversed(open(opt['input']).readlines()):
        if 'Optimization has failed' in line.strip():
            print("Psi4 failed optimization. Not extracting coordinates.")
            quit()
        if line.strip().startswith('Cleaning optimization helper files'):
            keeper = True
        if keeper:
            keeplines.append(line.strip())
        if line.strip().startswith('Saving final (previous) structure'):
            break
    keeplines.reverse()
    keeplines = keeplines[11:-2]
    print len(keeplines)

    if opt['type'].lower() == 'xyz':
        with open(opt['output'], 'w') as of:
            of.write('%d\n\n' % len(keeplines))
            for i in keeplines:
                of.write('%s\n' % i)
    else:
        atomLabels = ['N1 ','C2 ','N3 ','H4 ','H5 ','N6 ','C7 ','N8 ','C9 ','C10','C11','C12','C13','C14','N15','H16','H17','H18','H19','H20','H21','H22','H23']

        with open(opt['output'], 'w') as of:
            of.write('CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n')
            for i, atomLabel in zip(keeplines, atomLabels):
                trim = i[12:18]+' '+i[31:37]+' '+i[50:56]
                atom = atomLabel[0]
                count = atomLabel[-2:]

                of.write('ATOM     %s  %s GBI1X   1      %s  1.00  0.00      GBI1 %s\n' % (count, atomLabel, trim, atom))
            of.write('END\n')
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        help="Psi4 output filename (with path) from which to get final geometry.")
    parser.add_argument("-o", "--output",
                        help="Output filename. Example: ./gbi.xyz ")
    parser.add_argument("-t", "--type",
                        help="XYZ or PDB format. This script is for taut2 only!!")

    args = parser.parse_args()
    opt = vars(args)
    main(**opt)

