
# Purpose:
# Usage:    python convergedfep.py -d [/path/to/FEP-F/and/FEP-R] -n [num_blocks]
# Example:  python convergedfep.py -d ../../../ -n 10

# Where: /beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/1_dG/convergedfep.py

import numpy as np
import argparse
import re, os, glob

def numericalSort(value):
    """
    Parses some number. 5 would return ['5']. 5.4 would return ['5', '4'].
    For use in shorten_cat_fepout.
    """
    numbers = re.compile(r'(\d+)') # parses a given value
    parts = numbers.split(value)
    parts[1::2] = list(map(int, parts[1::2]))
    return parts


def tail( f, lines ):
    """
    https://stackoverflow.com/questions/136168/get-last-n-lines-of-a-file-with-python-similar-to-tail
    """

    if lines == 0:
        return []

    BUFSIZ = 1024
    f.seek(0, 2)
    remaining_bytes = f.tell()
    size = lines + 1
    block = -1
    data = []

    while size > 0 and remaining_bytes > 0:
        if remaining_bytes - BUFSIZ > 0:
            # Seek back one whole BUFSIZ
            f.seek(block * BUFSIZ, 2)
            # read BUFFER
            bunch = f.read(BUFSIZ)
        else:
            # file too small, start from beginning
            f.seek(0, 0)
            # only read what was not read
            bunch = f.read(remaining_bytes)

        bunch = bunch.decode('utf-8')
        data.insert(0, bunch)
        size -= bunch.count('\n')
        remaining_bytes -= BUFSIZ
        block -= 1

    return ''.join(data).splitlines()[-lines:]

def shorten_cat_fepout(fep_dir, outfile, num_blocks, cur_block):


    # calculate number of lines to keep
    num_lines_in_file = 2506    # for my 5 ns setup. change if not applicable! <===== check me
    lines_wanted = int((num_lines_in_file/num_blocks)*cur_block)
    print("Extracting {} lines".format(lines_wanted))

    # get list of all *.fepout file in this fep_dir
    fep_file = sorted(glob.glob(fep_dir+'/*.fepout'), key=numericalSort)

    # loop through all *.fepout files and write to the summary file
    with open(outfile, 'w') as output:
        print('Concatenating {}'.format(outfile))
        for fname in fep_file:
            f = open(fname, "rb")
            extracted_lines = tail(f, lines_wanted)
            for i in extracted_lines:
                output.write("{}\n".format(i))

    return outfile





def convergedfep(way, **kwargs):

    src = args.hdir.rstrip('//')
    hdir = src+'/FEP_{}/results'.format(way)
    outfile =  '{}_{}.fepout'.format('results', way)

    for i in range(1,args.numblocks+1):

        # directory setup
        print("Processing {} block(s) from the end...".format(i))
        newfile = "{:02}/".format(i)+outfile

        # file/directory checking
        if os.path.isfile(newfile):
            print("!!! WARNING: {} already exists".format(newfile))
        elif os.path.exists(hdir) == True:
            if not os.path.exists("{:02}".format(i)): os.mkdir("{:02}".format(i))
            shorten_cat_fepout(hdir, newfile, args.numblocks, i)
        else:
            print(os.getcwd())
            raise OSError("No such file or directory '{}'".format(hdir))



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--hdir",
                        help="Location with both FEP_F and FEP_R directories")
    parser.add_argument("-n", "--numblocks", type=int,
                        help="Separate each window's fepout into this many blocks. "
                             "This will correspond to the number of concatenated "
                             "output files for each of F and R directions.")

    args = parser.parse_args()
    opt = vars(args)

    convergedfep('F', **opt)
    convergedfep('R', **opt)
