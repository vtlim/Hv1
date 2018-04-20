
import sys
import numpy as np
import matplotlib.pyplot as plt

def vector_magnitude(forces):
    mags = []
    for i in forces:
        mags.append(np.linalg.norm(i))
    return np.asarray(mags)

def proc_pair_int(**kwargs):

    # read input file and extract time and force data
    tstep = []
    vforce = []
    eforce = []
    with open(opt['input']) as f:
        for line in f:
            if line.startswith("PAIR INTERACTION:"):
                l = line.split()
                tstep.append(int(l[3]))
                vforce.append(np.asarray([float(i) for i in l[5:8]]))
                eforce.append(np.asarray([float(i) for i in l[9:12]]))

    # calculate magnitude of each force vector
    vmags = vector_magnitude(vforce)
    emags = vector_magnitude(eforce)

    sys.path.insert(0,'/DFS-L/DATA/mobley/limvt/analysis')
    import plotXY
    plt.plot(plotXY.moving_average(vmags,1))
    plt.plot(plotXY.moving_average(emags,1))
#    plt.plot(vmags)
#    plt.plot(emags)
    plt.show()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        help="Name of the input NAMD log file.")
    parser.add_argument("--vdw", default=True,
                        help="True = process vdW component of force.")
    parser.add_argument("--elec", default=True,
                        help="True = process electrostatic force component.")

    args = parser.parse_args()
    opt = vars(args)
    proc_pair_int(**opt)

