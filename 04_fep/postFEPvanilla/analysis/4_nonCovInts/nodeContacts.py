
"""

DESCRIPTION: Process node-edge data files for contact interactions.

Input: Data files from VMD Tcl scripts by Eric Wong & Alfredo Freites.
Output: Pickle file with pandas dataframe with headers for the following:
  - node_i: nodes have their own indices, not corresponding to protein resid.
  - node_j
  - weight: number of contacts between two given nodes for a pair.
      I think every pair of a certain node_i and node_j should have same weight.
  - attribute: description of the edge based on node types in the pair.
      Attribute has the following options:
      HPHOB: between two nonpolar nodes
      COUL: between positive or negative nodes (4 combinations possible)
      HBOND: between a dipolar node and one of (dipolar, positive, negative)
      STER: the default, everything else
  - count: how many times (trajectory frames) the node pair had contact
  - average: weight/count

By: Victoria Lim, UCI

NOTES
* This is pretty intensive so would recommend interactive job request
  or Submitting the script as a job.
* Output is saved in a pickle file. Read this in via iPython notebook
  for analyzing and plotting data.

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import gc
import datetime


def readNewEdge(attFile, edFile):

    # Read the files in as pandas dataframes
    edatt = pd.read_csv(attFile,delimiter=' ',header=None)
    ednum = pd.read_csv(edFile,delimiter=' ',header=None,usecols=[2])

    # Combine the the attributes with the numbers, and add column titles
    edges = pd.concat([edatt,ednum],axis=1)
    edges.columns = ['node_i','node_j','attribute','weight']
    edges = edges.sort_values(['node_i','node_j'], ascending=[True,True])

    return edges


def mergeBySum(df1, df2):
    '''

    Joins data frames df1 and df2, keeping unique values as is
    and taking the sum for entries with matching node_i
    and node_j values.

    df1 (parent) should have columns: ['node_i','node_j','attribute','weight','count']
    df2 (new)    should have columns: ['node_i','node_j','attribute','weight']

    '''

    # combine with existing data
    bigdf = pd.concat((df1,df2))
    # sum the total number of counts of each node pair (pair is node_i to node_j)
    bigdf['count'] = bigdf['count'] + bigdf.groupby(['node_i', 'node_j'])['weight'].transform(np.size)-1
    # for node pairs listed multiple times, sum their weights
    filtered = bigdf.groupby(['node_i', 'node_j'],as_index=False)['weight'].sum()
    # all attributes and new pairs lost after summing so get them back
    filtered = filtered.merge(bigdf[['node_i', 'node_j','attribute','count']], how='inner', on=['node_i', 'node_j'])
    filtered = filtered.drop_duplicates(subset=['node_i', 'node_j','weight','attribute'])
    # fill in the NaN counts of new pairs to be 1.0
    filtered['count'].fillna(1.0, inplace=True)

    return filtered


def nodeContacts(nodeFile, edgePrefix, frame_a, frame_b, pickleOut):
    """
    """
    ### Read in node information.
    nodes = pd.read_csv(nodeFile,delimiter='\t',header=None,usecols=range(7))
    nodes.columns = ['index','resname','resid','location','type','code','nAtoms']
    nodes = nodes.set_index('index')
    #nodes.describe()


    ### Read in edges' information

    # Start the data frame with the first file's information
    edges = readNewEdge('{}_{}.edgeatt'.format(edgePrefix, frame_a),'{}_{}.edges'.format(edgePrefix, frame_a))
    edges['count'] = 1

    k = str(int(frame_a) + 1).zfill(len(frame_a)) # increment with leading zeroes
    while k <= frame_b:
        if int(k)%10==0:
            print(k)
            print(datetime.datetime.now())
        newEdges = readNewEdge('{}_{}.edgeatt'.format(edgePrefix, k) ,'{}_{}.edges'.format(edgePrefix, k))
        edges = mergeBySum(edges, newEdges) # merge new data into the existing frame
        k = str(int(k) + 1).zfill(len(k)) # increment with leading zeroes
        del newEdges
        junk = gc.collect()

    edges['average'] = edges['weight']/edges['count']
    #edges.describe()


    ### Save information to pickle
    pickle.dump((nodes, edges), open( pickleOut, "wb" ) )
    # read with nodes, edges = pickle.load( open( pickleIn, "rb" ))


### ------------------- Parser -------------------

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-n", required=True,
        help="Name of the node file. Suffix should be .nodes")

    parser.add_argument("-e", required=True,
        help="Prefix of the edge file before frame number. Don't include '_'")

    parser.add_argument("-a", required=True,
        help="First frame number of edge file to read. Include leading zeroes if present.")

    parser.add_argument("-b", required=True,
        help="Last frame number of edge file to read.")

    parser.add_argument("-o", required=True,
        help="Name of output file with pickled data.")

    args = parser.parse_args()
    nodeContacts(args.n, args.e, args.a, args.b, args.o)
