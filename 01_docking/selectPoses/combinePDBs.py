#!/usr/bin/env python
### Purpose: Combine PDBs, 2x.
     # (1) Combine protein and ligand, both uncentered. (use this to compute transformation matrix for centering)
     # (2) Combine full system and centered ligand.
### Example: 
# python combinePDBs.py -d /pub/limvt/hv1/01_docking/2_tautomer -s x_together --posesFromFile
# python combinePDBs.py -d /pub/limvt/hv1/01_docking/2_tautomer -s z_togetherCentered -l y_ligCentered

import os, sys, glob
import shutil
import argparse

# ------------------------- Functions ------------------------- #

def makedir(dirname):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    else:
        sys.exit("%s directory already exists. Exiting." % dirname)



# ------------------------- Script ---------------------------- #

def main(**kwargs):
    os.chdir(os.path.join(opt['dir'], "pdbqt/wt_ans"))
    subdir = os.path.join(opt['dir'], "pdbqt/wt_ans/%s" % opt['sub'] )
    makedir(subdir)

    ### get pose list from revclus file
    if opt['posesFromFile']:
        poselist = []
        f = open('configurations.dat','r')
        data = f.readlines()
        f.close()
    
        for line in data:
           poselist.append(line.split()[1])

    ### get pose list from list of centered ligands
    else:
        poselist = glob.glob(opt['dir']+'/pdbqt/wt_ans/'+opt['lig']+'/*.pdb')
        print(opt['dir']+'/'+opt['lig']+'/*.pdb')

    ### loop over the docked pose pdbs
    for pdb in poselist:
        print(pdb)

        ### Get the protein and lig IDs
        if opt['posesFromFile']:
            # Ex. hv1-fit-17480_bs1_docked_ligand_11.pdb
            config = pdb.split('-')[2].split('_')[0]
            pose = pdb.split('_')[4].split('.')[0]

        else:
            # Ex. hv1-config-19793_gbi-pose-12.pdb
            config = pdb.split('-')[2].split('_')[0]
            pose = pdb.split('-')[4].split('.')[0]

        ### Get PDB filenames for protein and for ligand
        if opt['posesFromFile']:
            crd1 = os.path.join(opt['dir'], "pdb/hv1-fit-%s.pdb" % config)
            crd2 = os.path.join(opt['dir'], "pdbqt/wt_ans/hv1-fit-%s/bs_1/hv1-fit-%s_bs1_docked_ligand_%s.pdb"\
 % (config, config, pose))
        else:
            crd1 = os.path.join(opt['dir'], "pdb/full-hv1-system/hv1_depolarized_frame_%s.pdb" % config)
            #crd2 = os.path.join(opt['dir'], "pdbqt/wt_ans/%s/%s" % (opt['lig'], pdb))
            crd2 = pdb

        ### Get data from both files
        f = open(crd1,'r')
        data1 = f.readlines()[:-1]  # skip 'END' line in protein file
        f.close()

        f = open(crd2,'r')
        if opt['posesFromFile']: data2 = f.readlines()[7:]    # skip 'REMARK' lines in ligand file
        else: data2 = f.readlines()[1:-1]    # skip header and footer in ligand file
        f.close()

        ### Write out both protein and ligand to one PDB file.
        outfile = os.path.join(opt['dir'], "pdbqt/wt_ans/%s/hv1-config-%s_gbi-pose-%s.pdb"\
 % (opt['sub'], config, pose))
        with open(outfile, 'w') as output:
            output.writelines(data1)
            output.writelines(data2)
            output.write("END\n")


# ---------------------- Parse Arguments ---------------------- #

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dir",
                        help="Location with both FEP_F and FEP_R directories")
    parser.add_argument("-s", "--sub",
                        help="Name of the subdir within dir (not full path)")
    parser.add_argument("-l", "--lig", default='',
                        help="Name of the subdir within dir (not full path)\
containing the centered ligands")
    parser.add_argument("-o", "--outfname", default='results',
                        help="Base name of saved plots")
    parser.add_argument("--posesFromFile", action="store_true", default=False,
                        help="Get list of poses from configurations.dat file,\
as opposed to globbing a directory.")

    args = parser.parse_args()
    opt = vars(args)
    main(**opt)

