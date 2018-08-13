
# Purpose:  Concatenate all the windows' fepout files, but only take last bit from each.
# Usage:    python convergedfep.py -d [/path/to/FEP-F/and/FEP-R] -n [num_blocks]
# Example:  python convergedfep.py -d ../../../ -n 10
# Where: /beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/1_dG/convergedfep.py
#
# Notes:
#   - If you plan to use a different directory other than FEP_*/results for the
#     .fepout files, MAKE SURE TO CHANGE THAT IN THIS SCRIPT. For example,
#     you may want to evaluate the extended jobs in FEP_*/results-10ns.
#


import numpy as np
import argparse
import re, os, glob
import subprocess

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

    Parameters
    ----------
    f : opened file object
        This is not the name of the file, but the object returned from Python open().
    lines : integer
        Number of lines desired from the tail end of the file.

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


def grab( f, lines, start ):
    """
    """

    if lines == 0:
        return []
    wholelines = f.readlines()
    finish = start+lines+1

    return wholelines[start:finish]


def shorten_cat_fepout(fep_dir, fepout_length, equil_length, outfile, num_blocks, cur_block, frombeg, isolated):

    def write_file(fobj, lines, binary=False):

        # write out the extracted lines
        if binary:
            for i in lines:
                fobj.write("{}\n".format(i))
        else:
            for i in lines:
                fobj.write("{}".format(i))



    # calculate how many lines of data contribute to ensemble average data
    data_length = fepout_length - equil_length

    # calculate number of lines to keep
    lines_wanted = int((fepout_length/num_blocks)*cur_block)
    print("Extracting {} lines".format(lines_wanted))

    # get list of all *.fepout file in this fep_dir
    fep_file = sorted(glob.glob(fep_dir+'/*.fepout'), key=numericalSort)

    # loop through all *.fepout files and write to the summary file
    with open(outfile, 'w') as output:
        print('Concatenating {}'.format(outfile))
        for fname in fep_file:

            if isolated:
                if frombeg:    # from beginning, by blocks
                    # todo
                    pass
                else:          # from end, by blocks
                    f = open(fname, "r")
                    # calculate number of lines to keep
                    lines_wanted = int(fepout_length/num_blocks)
                    start_line = int((num_blocks-cur_block)*lines_wanted+4) # 4 is approximate <===== MAKE MORE APPLICABLE
                    if lines_wanted < data_length:
                        output.write("#0 STEPS OF EQUILIBRATION AT LAMBDA 0.0 COMPLETED\n") # compatible with VMD parse FEP
                        output.write("#STARTING\n") # to be compatible with bar4fep.py
                    print("Extracting {} lines".format(lines_wanted))
                    extracted_lines = grab(f, lines_wanted, start_line)
                    write_file(output, extracted_lines, False)
                    if lines_wanted < fepout_length:
                        output.write("#Free energy change for lambda window [ 0 0 ] is 0.0 ; net change until now is 0.0\n") # compatible with VMD parse FEP


            else:
                # calculate number of lines to keep
                lines_wanted = int((fepout_length/num_blocks)*cur_block)
                print("Extracting {} lines".format(lines_wanted))

                if frombeg:    # from beginning, inclusive chunks
                    with open(fname) as f:
                        extracted_lines = [next(f) for x in range(lines_wanted)]
                    write_file(output, extracted_lines, False)
                    if lines_wanted < fepout_length:
                        output.write("#Free energy change for lambda window [ 0 0 ] is 0.0 ; net change until now is 0.0\n") # compatible with VMD parse FEP

                else:           # from end, inclusive chunks
                    f = open(fname, "rb")
                    if lines_wanted < data_length:
                        output.write("#0 STEPS OF EQUILIBRATION AT LAMBDA 0.0 COMPLETED\n") # compatible with VMD parse FEP
                        output.write("#STARTING\n") # to be compatible with bar4fep.py
                    extracted_lines = tail(f, lines_wanted)
                    write_file(output, extracted_lines, True)

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
            shorten_cat_fepout(hdir, args.feplength, args.equlength, newfile, args.numblocks, i, args.frombeg, args.isolated)
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
    parser.add_argument("-l", "--feplength", type=int, default=2506,
                        help="How many lines are in a single fepout file. "
                             "Default is set at 2506, for a 5 ns window with "
                             "an output frequency of 1000.")
    parser.add_argument("-e", "--equlength", type=int, default=504,
                        help="How many lines precede the data collection line. "
                             "This value should include the header saying "
                             "\"STARTING COLLECTION OF...\" Default is set at "
                             "504, for a 1 ns equil with output freq 1000.")
    parser.add_argument("--frombeg", default=False, action='store_true',
                        help="Do you want to evaluate convergence by taking "
                             "successive amounts of data from the end "
                             "(frombeg=False) or from the beginning (frombeg=True)?")
    parser.add_argument("--isolated", default=False, action='store_true',
                        help="If output file is broken into N chunks, do you "
                             "want to grab the last 1, 1-2, 1-3, etc. chunks "
                             "(isolated=False) or do you want to be able to "
                             "grab chunk 3 without grabbing 1 and 2? (isolated"
                             " = True)")

    args = parser.parse_args()
    opt = vars(args)

    convergedfep('F', **opt)
    convergedfep('R', **opt)
