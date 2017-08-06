#!/usr/bin/python
# Purpose: Take pdb from last FEP frame and process as new system.
#   This script removes disappearing atoms of hybrid sys, renames mutated
#   residue to the final residue name, and removes the 'B' suffix on
#   appearing atom names.
# usage:   python script.py in.pdb out.pdb mutation finalres
# example: python script.py in.pdb out.pdb D2E GLU
# todo: make it more efficient than looping over indiv lines?

import sys

def editFEPpdb(inputf, outf, mutcode, fincode):

    deleted = 0
    renamed = 0
    with open(inputf,'r') as f, open(outf,'w') as g:
        for line in f:
            if mutcode in line:
                lcopy = line.split()
                oldatom = lcopy[2]
                # ignore starting residue represented by NONbackbone
                # atom name ending with A
                if len(oldatom) > 2 and oldatom[-1] == 'A':
                    deleted += 1
                    continue
                elif len(oldatom) > 2 and oldatom[-1] == 'B':
                    newatom = oldatom[:-1]+' '
                    newline = line.replace(mutcode, fincode)
                    newline = newline.replace(oldatom, newatom)
                    renamed += 1
                    g.write(newline)
                else:
                    newline = line.replace(mutcode, fincode)
                    renamed += 1
                    g.write(newline)
            else:
                g.write(line)

    print("Number of removed lines: %d" % (deleted))
    print("Number of altered lines: %d" % (renamed))



if __name__ == "__main__":
    editFEPpdb(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
