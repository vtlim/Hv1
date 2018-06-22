
"""
Analysis functions for nodeAnalysis.ipynb.
Victoria T. Lim
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import seaborn as sns



def getEdgePair(node1, node2, edges):
    """
    Get edge information that involves the specified two nodes.
    This function uses node indices as parameters; to pass in
    protein residue indices, use the function getContactPair.

    Parameters
    ----------
    node1 : int
        Index of node_i. This is the node index not protein index.
    node2 : int
        Index of node_j. Should be larger than node_i.
    edges : pandas dataframe
        Dataframe in which to search for interaction.

    Returns
    -------
    pandas dataframe

    Notes
    -----
    Can't find the interaction you're looking for? Might need to
    add/subtract an offset if you used the diffEdges function.
    The offset value (with sign and magnitude) should have been
    printed out: "Shifting node indices by..."

    """
    return edges[(edges.node_i == node1) & (edges.node_j == node2)]


def getContactPair(res1, res2, nodes, edges):
    """
    Get edge information that involves the two specified protein residues.
    This function is similar to getEdgePair but takes into protein indices.

    Parameters
    ----------
    res1 : int
        protein residue number to use for node_i
    res2 : int
        protein resiude number to use for node_j. Should be greater than res1.
    nodes : pandas dataframe
        Dataframe in which to translate protein index to node index (indices).
    edges : pandas dataframe
        Dataframe in which to search for interaction.

    Returns
    -------
    pandas dataframe

    """
    df1 = getResidInfo(res1,nodes,resExcludes=['WAT'])
    df2 = getResidInfo(res2,nodes,resExcludes=['WAT'])
    indexList1 = df1.index.tolist()
    indexList2 = df2.index.tolist()
    print(df1, '\n\n', df2)

    return edges[(edges.node_i.isin(indexList1)) & (edges.node_j.isin(indexList2))]


def findInEdges(nodeNum, edges, att=None):
    """
    Find the specified node index in either node_i or node_j columns of input edges df.

    Parameters
    ----------
    nodeNum : int
        This is the node index not protein index.
    edges : pandas dataframe
        Dataframe in which to search for node.

    Returns
    -------
    pandas dataframe

    """
    edges = edges.copy()
    if att is not None:
        return edges[(edges.attribute == att) & ((edges.node_i == nodeNum) | (edges.node_j == nodeNum))]
    else:
        return edges[(edges.node_i == nodeNum) | (edges.node_j == nodeNum)]


def getResidInfo(resid, nodes, resExcludes=[]):
    """
    Get the node information for specified protein residue index.

    Parameters
    ----------
    resid : int
        protein residue number
    nodes : pandas dataframe
        Dataframe in which to translate protein index to node index (indices).
    resExcludes: list
        List containing strings for residues to ignore. E.g., ['WAT']

    Returns
    -------
    pandas dataframe

    """
    nodes_id = nodes.loc[nodes['resid'] == resid]
    nodes_id = nodes_id[~nodes_id.resname.isin(resExcludes)]

    return nodes_id


def idxToResid(idx, nodes, idOnly=False):
    """
    This function takes in some node index and generates
    a string code of one-letter residue name and integer of residue number.

    Parameters
    ----------
    idx : int
        integer index of the pandas dataframe
    nodes : pandas dataframe
        pandas dataframe of which to search
    idOnly : Boolean
        True to return numpy.int64 of residue number
        False to return code with resname abbrev
        ex., True returns 150; False returns 'F150:sc'

    Returns
    -------
    string or numpy.int64 value of residue (based on idOnly parameter)

    """
    aa_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'HSD': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
     'GBI1':'GBI1', 'GBI2':'GBI2', 'WAT':'WAT'}

    #old way not conducive to taking diff of dataframes
    #entry = nodes.iloc[idx-1]
    entry = nodes.loc[nodes.index == idx] # iloc gives series, loc gives dataframe
    entry = entry.T.squeeze() # convert from dataframe to series
    resname = entry['resname']
    if resname in ['GBI1','GBI2']:
        code = aa_dict[resname]+':'+entry['code']
    elif resname == 'WAT':
        code = aa_dict[resname]+str(entry['resid'])
    else: # if not GBI or WAT, must be Hv1
        if idOnly: code = entry['resid']
        else: code = aa_dict[resname]+str(entry['resid'])+':'+entry['location']

    return code


def trimEdgesTM(nodes, edges):
    """
    Process edges to remove edges that don't include transmembrane (TM)
    protein residue. In other words, only keep edges that involve at
    least one TM residue. The TM residue may contact a non-TM residue,
    water, or another TM residue.

    Parameters
    ----------
    nodes : pandas dataframe
        Pandas dataframe with information on nodes (residues)
    edges : pandas dataframe
        Pandas dataframe with information on edges (contacts)

    Returns
    -------
    pandas dataframe

    """

    # define which residues are in TM region
    seg1 = list(range(99,126))
    seg2 = list(range(134,161))
    seg3 = list(range(168,192))
    seg4 = list(range(198,221))
    segtm = seg1+seg2+seg3+seg4

    # get node indices for residues in TM region
    protein_nodes = nodes[(nodes['resid'].isin(segtm)) & (nodes['resname'] != 'WAT')]
    prot_node_ids = protein_nodes.index.tolist()

    # keep edges with at least one node that is a TM residue
    return edges[ (edges['node_i'].isin(prot_node_ids)) | (edges['node_j'].isin(prot_node_ids)) ]


def prioritize(edges, rawNum):
    """
    TODO
    Pull out N strongest interactions, or pivot the table.
    Not meant for user; implemented in protLigInts and selectionInts functions.
    This function handles cases of whether dataframe has edges or difference of edges.

    Parameters
    ----------
    edges : pandas dataframe
    rawNum : integer

    Returns
    -------
    pandas dataframe

    """
    edges = edges.copy()
    try: # Pull out the N strongest interactions (should be no negative values)
        edges = edges.sort_values('average',ascending=False).head(rawNum)
    except KeyError: # if no 'average' column then this is one df minus another so there are negatives
        # sort by magnitude to get + and - changes
        tempinds = edges.avg_subt.abs().sort_values(ascending=False).head(rawNum).index
        edges = edges.loc[tempinds]
    return edges


def pivot(edges, data=""):
    """
    TODO
    Pull out N strongest interactions, or pivot the table.
    Not meant for user; implemented in protLigInts and selectionInts functions.
    This function handles cases of whether dataframe has edges or difference of edges.

    Parameters
    ----------
    edges : pandas dataframe
    rawNum : integer

    Returns
    -------
    pandas dataframe

    """
    edges = edges.copy()
    if data=="edgetype":
        edges = edges.pivot(index='node_i',columns='node_j', values='edgetype')
    else:
        try:
            edges = edges.pivot(index='node_i',columns='node_j', values='average')
        except KeyError:
            edges = edges.pivot(index='node_i',columns='node_j', values='avg_subt')
    edges = edges.dropna(axis=1,how='all') # drop columns with all nan's

    return edges


def protLigInts(nodes, edges, rawNum=250, dry=1):
    """
    Take in a set of nodes and edges and identify the N strongest interactions.
    This function disregards:
      (1) interactions between waters (there can be protein-water interaction),
      (2) interactions between adjacent residues (e.g., residue F149 and F150), and
      (3) interactions within the same residue (e.g., backbone and sidechain of F150).

    Parameters
    ----------
    nodes : pandas dataframe
        Pandas dataframe with information on nodes (residues)
    edges : pandas dataframe
        Pandas dataframe with information on edges (contacts)
    rawNum : integer
        How many interactions to use before further processing.
        Further processing = remove adjacent & intra-residue interactions.
    dry : integer
        2 means no waters at all even to protein/ligand
        1 means no water-water interactions
        0 means allow waters (NOT YET implemented)

    Returns
    -------
    pandas PIVOTED dataframe with reduced and filtered interactions, formatted as:
     > node_i as index column
     > node_j as different columns
     > average interaction strength in cell intersecting node_i and node_j

    """
    edges = edges.copy()

    # Get all indices of nodes that are not water
    watless_idx = nodes.index[nodes['resname'] != 'WAT'].tolist()
    if dry==1:
        # Filter interactions with at least one non-water (remove wat-wat interactions)
        watless_edges = edges.loc[edges['node_i'].isin(watless_idx) | edges['node_j'].isin(watless_idx)]
    elif dry==2:
        # Filter interactions with no waters whatsoever
        watless_edges = edges.loc[edges['node_i'].isin(watless_idx) & edges['node_j'].isin(watless_idx)]


    # Pull out the N strongest interactions
    watless_edges = prioritize(watless_edges,rawNum=rawNum)
    if watless_edges is None: return

    # Make temp copy to compare protein resIDs to filter out those in same/adj resid
    temp = watless_edges.copy()
    temp['node_i'] = temp['node_i'].apply(idxToResid,args=(nodes,True))
    temp['node_j'] = temp['node_j'].apply(idxToResid,args=(nodes,True))
    # convert the non-protein residues with no ID for temp integer
    temp['node_i'].replace('GBI\w', -500, regex=True,inplace=True)
    temp['node_j'].replace('GBI\w', -500, regex=True,inplace=True)
    temp['node_i'].replace('WAT\w', -400, regex=True,inplace=True)
    temp['node_j'].replace('WAT\w', -400, regex=True,inplace=True)

    # drop node interactions in same resid or adjacent
    dropinds = temp.index[((temp['node_i']-temp['node_j']).abs() <= 1) == True].tolist()
    watless_edges.drop(dropinds, inplace=True)

    return watless_edges


def selectionInts(nodes, edges, indices, rawNum=50, dry=True):
    """
    Parameters
    ----------
    nodes : pandas dataframe
        Pandas dataframe with information on nodes (residues)
    edges : pandas dataframe
        Pandas dataframe with information on edges (contacts)
    indices : list of integers
        List of node indices of selection. Two examples:
        1. gidx_1 = nodes_1.index[nodes_1['resname'] == 'GBI1'].tolist()
        2. selNodes = getResidInfo(211, nodes_2, resExcludes=['WAT'])
           selInds = selNodes.index.tolist()
    rawNum : int
        How many interactions to use before further processing.
        Further processing = remove adjacent & intra-residue interactions.
    dry : Boolean
        True to ignore any water-interactions of given selection

    Returns
    -------
    sel_edges - pandas PIVOTED dataframe with interactions for given selection
                new format:
                 > node_i as index column
                 > node_j as different columns
                 > average interaction strength in cell intersecting node_i and node_j


    Examples of selecting indices
    -----------------------------
    > gidx_1 = nodes_1.index[nodes_1['resname'] == 'GBI1'].tolist()

    > selNodes = getResidInfo(211, nodes_2, resExcludes=['WAT'])
    > selInds = selNodes.index.tolist()

    """
    sel_edges = edges.copy()

    if dry:
        watidx = nodes.index[nodes['resname'] == 'WAT'].tolist() # get water indices
        sel_edges = sel_edges[(~sel_edges['node_i'].isin(watidx)) & (~sel_edges['node_j'].isin(watidx))]

    # Get all the edge interactions that relate to selection
    sel_edges = sel_edges.loc[sel_edges['node_i'].isin(indices) | sel_edges['node_j'].isin(indices)]

    # Pull out the N strongest interactions
    sel_edges = prioritize(sel_edges,rawNum=rawNum)
    if sel_edges is None: return

    # put all the GBI nodes in the i spot
    sel_edges["node_i"], sel_edges["node_j"] = np.where(sel_edges['node_j'].isin(indices),
        [sel_edges["node_j"], sel_edges["node_i"]], [sel_edges["node_i"], sel_edges["node_j"]])
    sel_edges = sel_edges[~sel_edges['node_j'].isin(indices)] # remove self-interactions

    return sel_edges


def plotHeatInts(nodes,edges,minHeat=0,maxHeat=20,colors=None,size=(20,20),seltitle="",pivoted=False):
    """

    Parameters
    ----------
    nodes : pandas dataframe
    edges : pandas dataframe
        Pivoted pandas dataframe (node_i indices as header, node_j indices as left column.)
    minHeat : integer
        Minimum data point in heat color bar.
        Should be zero unless there's a special case.
    maxHeat : integer
        Maximum data point in heat color bar.
        May want to manually adjust if max edge data > default maxHeat.
    colors : string
        String code referring to one of Python's color maps.
        https://matplotlib.org/examples/color/colormaps_reference.html
        Use a parameter of colors="edges" to color by edge type instead of strength.
    size : tuple
        Tuple of length 2 for (width, length) in matplotlib
    seltitle : string
        Title to list at the top of the plot
    pivoted : Boolean
        Whether or not the input edges is already pivoted (such as from s

    """
    def offsetHeatGrid():
        # offset the y-grid to match the label WITHOUT offsetting ticklabels
        yticks_old = ax.get_yticks()
        if len(yticks_old) > 1:
            yticks_offset = (yticks_old[1]-yticks_old[0])/2
            yticks = [(tick-yticks_offset) for tick in ax.get_yticks()]
            ax.set_yticks(yticks) # grid will use these new tick placements
            ax.set_yticks(yticks_old,minor=True)
            ax.set_yticklabels(ylabels,minor=True) # put labels back in old placements
            ax.set_yticklabels([]) # turn off labels at new tick placements

        # offset the x-grid to match the label WITHOUT offsetting ticklabels
        xticks_old = ax.get_xticks()
        if len(xticks_old) > 1:
            xticks_offset = (xticks_old[1]-xticks_old[0])/2
            xticks = [(tick-xticks_offset) for tick in ax.get_xticks()]
            ax.set_xticks(xticks) # grid will use these new tick placements
            ax.set_xticks(xticks_old,minor=True)
            ax.set_xticklabels(xlabels,minor=True) # put labels back in old placements
            ax.set_xticklabels([]) # turn off labels at new tick placements

    def label_edge(row):
        # https://stackoverflow.com/questions/26886653/pandas-create-new-column-based-on-values-from-other-columns
        if row['attribute'] == "HPHOB": # two nonpolar nodes
            return 1
        if row['attribute'] == "COUL":  # two charged nodes
            return 2
        if row['attribute'] == "HBOND": # dipolar and charged nodes
            return 3
        if row['attribute'] == "STER":  # everything else
            return 4
        return -1


    if (colors=="edgetype" and pivoted==True):
        print("Cannot color by edgetype if you have pass in a pivoted edge plot.")
        return

    if pivoted:
        plotInput = edges
    elif colors=="edgetype":
        # reassign the strength values in the edges dataframe
        plotInput = edges.copy()
        plotInput['edgetype'] = plotInput.apply (lambda row: label_edge(row),axis=1)
        plotInput = pivot(plotInput, data='edgetype')
    else:
        plotInput = pivot(edges)

    plotNodes = nodes

    # generate plot labels based on residue name and residue number
    ylabels = [idxToResid(i, plotNodes) for i in list(plotInput)] # get node_j's, convert idxToResid
    xlabels = [idxToResid(i, plotNodes) for i in list(plotInput.index.values)] # get node_i's, convert idxToResid

    # plot the data
    plt.clf()
    plt.subplots(figsize=size)
    sns.set(font_scale=2.1)
    if colors=='edgetype':
        colors="tab10"
        vmin=0
        vmax=4
    ax = sns.heatmap(plotInput.T,annot=True,yticklabels=ylabels,xticklabels=xlabels,
                     cmap=colors,vmin=minHeat, vmax=maxHeat)

    offsetHeatGrid()

    plt.grid()
    plt.ylabel('')
    plt.xlabel('')
    plt.yticks(rotation=0)
    plt.title('\"Strongest\" interactions of {}'.format(seltitle))
    plt.show()
    print('interaction range is from {} to {}; verify if this is appropriate'.format(minHeat,maxHeat))


def diffEdges(nodes_x,nodes_y,edges_x,edges_y):
    """
    Take one set of edges and subtract another. This can identify changes in contacts
    between different systems, such as before and after mutation.

    USE CASES:
       [1] taut1 and taut2 with SAME protein system but differs in 2GBI and maybe in waters
       [2] tautx before and after mutation, all else the same (same 2GBI and waters)

    """
    nodes_x = nodes_x.copy()
    nodes_y = nodes_y.copy()
    edges_x = edges_x.copy()
    edges_y = edges_y.copy()

    # take union of both dataframes wrt to nodes_y, and ...
    df_1 = pd.merge(nodes_y, nodes_x, how='outer', indicator=True)
    # do again to get orig indices of nodes_x, ...
    df_2 = pd.merge(nodes_x, nodes_y, how='outer', indicator=True)
    # in order to find rows not in one or other
    rows_in_df1_not_in_df2 = df_1[df_1['_merge']=='left_only'][nodes_y.columns]
    rows_in_df2_not_in_df1 = df_2[df_2['_merge']=='left_only'][nodes_x.columns]
    # convert those rows to list of indices
    nodes_in_df1_not_in_df2 = np.asarray(rows_in_df1_not_in_df2.index.tolist())+1
    nodes_in_df2_not_in_df1 = np.asarray(rows_in_df2_not_in_df1.index.tolist())+1
    print("nodes in 1st, not in 2nd: ",nodes_in_df1_not_in_df2)
    print("nodes in 2nd, not in 1st: ",nodes_in_df2_not_in_df1)

    # remove the those indices from each set of edges respectively
    noMut_edges_y = edges_y.loc[~edges_y['node_i'].isin(nodes_in_df1_not_in_df2)]
    noMut_edges_y = noMut_edges_y.loc[~noMut_edges_y['node_j'].isin(nodes_in_df1_not_in_df2)]
    noMut_edges_x = edges_x.loc[~edges_x['node_i'].isin(nodes_in_df2_not_in_df1)]
    noMut_edges_x = noMut_edges_x.loc[~noMut_edges_x['node_j'].isin(nodes_in_df2_not_in_df1)]


    # scale edges_x to match indices of nodes_y
    mut1_start = nodes_in_df1_not_in_df2[0]
    mut2_start = nodes_in_df2_not_in_df1[0]
    if not mut1_start == mut2_start:
        print("ERROR: start node value of mutation(?) do not match between nodes_x and nodes_y")
        sys.exit()
    mut1_end = nodes_in_df1_not_in_df2[-1]
    mut2_end = nodes_in_df2_not_in_df1[-1]
    offset = mut1_end - mut2_end # NOT valid if there are 2 sets of changes (e.g., mutation && taut)
    if np.absolute(offset) > 10:
        print("ERROR: nodes are offset by {}. Check use cases, tail of nodes, and nodes in one but not the other.".format(offset))
        sys.exit()
    mask_i = (noMut_edges_x['node_i'] > mut2_end)
    mask_j = (noMut_edges_x['node_j'] > mut2_end)
    to_change_i = noMut_edges_x[mask_i]
    to_change_j = noMut_edges_x[mask_j]

    print("Shifting node indices by {} for {} rows".format(offset,len(to_change_i.index)+len(to_change_j.index)))
    to_change_i['node_i'] += offset
    to_change_j['node_j'] += offset
    new_nodes_i = to_change_i.pop('node_i')
    new_nodes_j = to_change_j.pop('node_j')
    noMut_edges_x.loc[noMut_edges_x.index.isin(to_change_i.index), 'node_i'] = new_nodes_i
    noMut_edges_x.loc[noMut_edges_x.index.isin(to_change_j.index), 'node_j'] = new_nodes_j

    # subtract the set of edges based on matching nodes (x-y)
    diff_edges_x = noMut_edges_x.merge(noMut_edges_y,on=['node_i','node_j','attribute'],how='inner')
    diff_edges_x['avg_subt'] = diff_edges_x['average_x']-diff_edges_x['average_y']

    # get nodes in common between both dataframes but keep nodes_y indices, https://tinyurl.com/y9qwemhv
    nodes_unified = nodes_y.reset_index().merge(nodes_x, how='inner', on=['resname','resid','location','type','code','nAtoms']).set_index('index')

    return nodes_unified, diff_edges_x, mut1_start, offset


def diffdiffEdges(nodes_x,nodes_y,edges_x,edges_y):
    """
    Similar to diffEdges function but the input edges here have already been through diffEdges.
    Processes the input edges to remove *_x and *_y columns and rename "avg_subt" to "average" column.
    This function does call the diffEdges function.

    Parameters
    ----------

    Returns
    -------
    """
    de_x = edges_x.copy()
    de_y = edges_y.copy()

    ### step 1: dry_diff_edges_1
    # process edges df to remove weight_x, count_x, average_x (and *_y); rename "avg_subt"
    de_x = de_x[['node_i','node_j','attribute','avg_subt']]
    de_x.rename(columns={'avg_subt':'average'}, inplace=True)

    ### step 2: dry_diff_edges_0
    # process edges df to remove weight_x, count_x, average_x (and *_y); rename "avg_subt"
    de_y = de_y[['node_i','node_j','attribute','avg_subt']]
    de_y.rename(columns={'avg_subt':'average'}, inplace=True)

    ### step 3: dry_diff_edges_1 - dry_diff_edges_0
    # the "nodes in 1st, not in 2nd" results make less sense here since diff of diff
    diff_nodes_10, diff_edges_10, mutstart_10, offset_10 = diffEdges(nodes_x,nodes_y,de_x,de_y)

    return diff_nodes_10, diff_edges_10, mutstart_10, offset_10


def similarizeTwoEdges(edges_x, edges_y, mutstart=None, offset=None, ignoreI=True):
    """
    Take in two similar dataframes of edges and find common interactions.
    In other words, extract common node_i's and node_j's of each.

    TODO: MORE IMPORTANT: rewrite this function. Originally was written for
        two pivoted dataframes that were passed in, but since pivoting is now
        done right before plotting (to keep edgetype information), this script
        should be updated. After that, you should be able to remove the "pivoted"
        variable in the plotHeatInts function.
    TODO: make segregation of i and j better. right now, checking ignoreI only in find common
        section, but if ignoreI is true, don't need to do ANY of the row stuff.

    Parameters
    ----------
    edges_x
    edges_y
    mutstart | int | integer value of [where] the two node sets start to become offset if at all
                     This should come directly from the diffEdges function, else determine manually.
                     If two sets of nodes are identical, leave as None.
    offset   | int | integer value of [how much] the two node are offset starting after mutstart
                     This should come directly from the diffEdges function, else determine manually.
                     If two sets of nodes are identical, leave as None.
                     For example, if GBI starts at node 11231 for edges_x and at node 11230 for edges_y,
                     then the offset is (x-y) for an integer value of 1.
    ignoreI | Bool | True means don't worry about making the node_i's uniform, just the node_j's

    Returns
    -------
    edges_x
    edges_y

    """
    edges_x = pivot(edges_x)
    edges_y = pivot(edges_y)

    # get list of ROWS (node_i's) of each df. copy to not chg orig idx.
    is_from_x = np.copy(np.asarray(edges_x.index.values))
    is_from_y = np.copy(np.asarray((edges_y.index.values)))

    # get list of COLUMNS (node_j's) of each df. copy to not chg orig idx.
    js_from_x = np.copy(np.asarray(edges_x.columns.values))
    js_from_y = np.copy(np.asarray(edges_y.columns.values))

    # if either is defined, then both must be defined
    if (mutstart or offset) and not (mutstart and offset):
        print("ERROR: Define both values for mutstart and offset.")
        return None, None

    # scale all node_j values of edges_y by offset if after mutstart.
    if mutstart is not None:
        for n, i in enumerate(is_from_y):
            if i>mutstart: is_from_y[n] = i+offset
        for n, j in enumerate(js_from_y):
            if j>mutstart: js_from_y[n] = j+offset

    # take node_j's in common for both set of edges
    common_is = np.intersect1d(is_from_x, is_from_y)
    common_js = np.intersect1d(js_from_x, js_from_y)
    print("node_i's in common: ", common_is)
    print("node_j's in common: ", common_js)

    # filter dataframes based on ROW indices in common
    if not ignoreI:
        edges_x = edges_x[edges_x.index.isin(common_is)]
        edges_y = edges_y[edges_y.index.isin(common_is)]

    # filter dataframes based on COLUMN indices in common
    edges_x = edges_x[common_js]
    edges_y = edges_y[common_js]

    return edges_x, edges_y


def condenseWaters(nodes, edges):
    """
    For each non-water node, sum up all the interaction weights involving water,
      and summarize this value into a single new edge.
    If used, call this function before calling protLigInts or selectionInts.

    Brief explanation
    -----------------
    Let's say you start with:
      > 11369 total original edges
      > 7743 edges have at least one water
      > 445 nodes that are not water
    Then add a new edge interaction for every not-water node.
    The average/avg_subt value is the sum of all waters.
      > 11814 (11369 + 445) total edges
    Then remove all the original edges that involve at least one water.
    This does not affected newly added edges bc watidx doesn't have index of new sumWat node.
      > 4071 (11814 - 7743) final edges

    """

    nodes = nodes.copy()
    edges = edges.copy()

    # separate indices to "water" and "not water:
    watidx = nodes.index[nodes['resname'] == 'WAT'].tolist()
    notidx = nodes.index[nodes['resname'] != 'WAT'].tolist()
    print("The water indices range from {} to {}".format(watidx[0],watidx[-1]))

    # node_j of the summed water is (value of last current node)+1
    newWatNode = nodes.index[-1]+1
    print("The new index of the summed waters is",newWatNode)

    # loop over each non-water node and SUM all water interactions into new edge row
    for n in notidx:
        # take subset of edges which involve this residue only
        subEdges = findInEdges(n,edges)
        try:
            iwat_sum = subEdges.loc[subEdges['node_i'].isin(watidx), 'average'].sum()
            jwat_sum = subEdges.loc[subEdges['node_j'].isin(watidx), 'average'].sum()
        except KeyError:
            iwat_sum = subEdges.loc[subEdges['node_i'].isin(watidx), 'avg_subt'].sum()
            jwat_sum = subEdges.loc[subEdges['node_j'].isin(watidx), 'avg_subt'].sum()
        totwat_sum = iwat_sum + jwat_sum
        #print(n, totwat_sum) # to troubleshoot if plots look odd

        # create new row in the main edges dataframe
        # node_i, node_j, weight_x, attribute, count_x, average_x, weight_y, count_y, average_y, avg_subt
        edges.loc[edges.shape[0]] = [n,newWatNode,'','HBOND','','','','','',totwat_sum]


    # add new node for the sum_wat with node ID of 0 (not neg. so regex can read in protLigInts)
    # resname, resid, location, type, code, nAtoms
    nodes.loc[newWatNode] = ['WAT','0','wt','DIP','OH2','1']

    # remove all the edges involving water except for new row
    edges = edges.loc[~edges['node_i'].isin(watidx) | ~edges['node_j'].isin(watidx)]

    return nodes, edges

def plotBarWaters(nodes, edges, rawNum=20, size=(40,10),pivoted=False):
    """
    Plot the interaction strength of the summed waters by residue.
    (Call this after condenseWaters function.)
    The edges dataframe are pivoted (if not already) and plotted as a bar chart.

    Parameters
    ----------
    nodes
    edges
    rawNum
    size
    pivoted

    Returns
    -------

    """

    if not pivoted:
        sumwatidx = nodes[nodes.resname=='WAT'].index.tolist()[-1]
        waterEdges = findInEdges(sumwatidx, edges)
        # place the node for summaryWater in the node_i spot
        waterEdges["node_i"], waterEdges["node_j"] = np.where(waterEdges['node_j'] == sumwatidx,
            [waterEdges["node_j"], waterEdges["node_i"]], [waterEdges["node_i"], waterEdges["node_j"]])
        waterEdges = prioritize(waterEdges,rawNum=rawNum)
        # after sort by value, sort the x-axis back by index again
        waterEdges = waterEdges.sort_values('node_j')
        # pivot table for plotting
        pivoted_waterEdges = pivot(waterEdges)
    else:
        pivoted_waterEdges = edges
    # plot
    xlabels = [idxToResid(i, nodes) for i in list(pivoted_waterEdges.columns.values)] # node_j's
    fig, ax = plt.subplots()
    ax = pivoted_waterEdges.T.plot(kind='bar',figsize=size,ax=ax,legend=False) # transpose to get columns as x-axis
    ax.set_xticklabels(xlabels)
    plt.show()

    if not pivoted:
        return waterEdges
    else:
        return pivoted_waterEdges


