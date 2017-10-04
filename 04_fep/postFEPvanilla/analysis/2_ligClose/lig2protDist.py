
# usage: mpython lig2protDist.py t1_v178a

import mdtraj as md
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys

key=sys.argv[1]

print('loading in {} trajectory...'.format(key))
t=md.load(key+'.dcd',top=key+'.psf',stride=100)

print('assigning residue groups...')
top = t.topology
gbiAtoms = top.select("resname GBI1")
group0 = [top.atom(gbiAtoms[0]).residue.index]

### mdtraj starts reading residue 0 as the first one (Phe 88)
# so to get mdtraj values, take (protein resid)-88
group1 = list(range(11,38)) # 99-126
group2 = list(range(46,73)) # 134-161
group3 = list(range(80,104)) # 168-192
group4 = list(range(110,133)) # 198-221


print('calculating contacts from s1 helix...')
pairs1 = list(itertools.product(group0, group1))
dists1, inds1 = md.compute_contacts(t,pairs1) # inds1 is same as pairs1 here
print('calculating contacts from s2 helix...')
pairs2 = list(itertools.product(group0, group2))
dists2, inds2 = md.compute_contacts(t,pairs2)
print('calculating contacts from s3 helix...')
pairs3 = list(itertools.product(group0, group3))
dists3, inds3 = md.compute_contacts(t,pairs3)
print('calculating contacts from s4 helix...')
pairs4 = list(itertools.product(group0, group4))
dists4, inds4 = md.compute_contacts(t,pairs4)

### take relative to reference coordinates
print('doing the same with reference coordinates...')
u=md.load('t1_begin.pdb')
group0a = [u.topology.atom(u.topology.select("resname GBI1")[0]).residue.index]
# assign pairs
pairs1a = list(itertools.product(group0a, group1))
pairs2a = list(itertools.product(group0a, group2))
pairs3a = list(itertools.product(group0a, group3))
pairs4a = list(itertools.product(group0a, group4))
# compute distances
dists1a, inds1a = md.compute_contacts(u,pairs1a)
dists2a, inds2a = md.compute_contacts(u,pairs2a)
dists3a, inds3a = md.compute_contacts(u,pairs3a)
dists4a, inds4a = md.compute_contacts(u,pairs4a)
# take relative difference
rel1 = dists1-dists1a
rel2 = dists2-dists2a
rel3 = dists3-dists3a
rel4 = dists4-dists4a


print('plotting original distances...')
plt.clf()
fig = plt.figure(figsize=(12,24))

plt.subplot(4,1,1)
plt.imshow(dists1.T,cmap='jet_r',aspect='auto',interpolation='none', vmin=0.0, vmax=2.30)
plt.ylabel('residue in S1')
ax = plt.gca();
ax.set_yticks(np.arange(0, 27, 1));
ax.set_yticklabels(np.arange(99, 126, 1));

plt.subplot(4,1,2)
plt.imshow(dists2.T,cmap='jet_r',aspect='auto',interpolation='none', vmin=0.0, vmax=2.30)
plt.ylabel('residue in S2')
ax = plt.gca();
ax.set_yticks(np.arange(0, 27, 1));
ax.set_yticklabels(np.arange(134, 161, 1));

plt.subplot(4,1,3)
plt.imshow(dists3.T,cmap='jet_r',aspect='auto',interpolation='none', vmin=0.0, vmax=2.30)
plt.ylabel('residue in S3')
ax = plt.gca();
ax.set_yticks(np.arange(0, 24, 1));
ax.set_yticklabels(np.arange(168, 192, 1));

plt.subplot(4,1,4)
plt.imshow(dists4.T,cmap='jet_r',aspect='auto',interpolation='none', vmin=0.0, vmax=2.30)
plt.ylabel('residue in S4')
ax = plt.gca();
ax.set_yticks(np.arange(0, 23, 1));
ax.set_yticklabels(np.arange(198, 221, 1));

plt.xlabel('every 100th frame')
plt.xticks(np.arange(0, 50, 5))
colorbar = plt.colorbar(cax = fig.add_axes([0.95, 0.12, 0.03, 0.76])) # left bot wid height
plt.savefig(key+'_byChain.eps',bbox_inches='tight')



print('plotting relative distances...')
plt.clf()
fig = plt.figure(figsize=(12,24))

plt.subplot(4,1,1)
plt.imshow(rel1.T,cmap='seismic_r',aspect='auto',interpolation='none', vmin=-0.5, vmax=0.5)
plt.ylabel('residue in S1')
ax = plt.gca();
ax.set_yticks(np.arange(0, 27, 1));
ax.set_yticklabels(np.arange(99, 126, 1));

plt.subplot(4,1,2)
plt.imshow(rel2.T,cmap='seismic_r',aspect='auto',interpolation='none', vmin=-0.5, vmax=0.5)
plt.ylabel('residue in S2')
ax = plt.gca();
ax.set_yticks(np.arange(0, 27, 1));
ax.set_yticklabels(np.arange(134, 161, 1));

plt.subplot(4,1,3)
plt.imshow(rel3.T,cmap='seismic_r',aspect='auto',interpolation='none', vmin=-0.5, vmax=0.5)
plt.ylabel('residue in S3')
ax = plt.gca();
ax.set_yticks(np.arange(0, 24, 1));
ax.set_yticklabels(np.arange(168, 192, 1));

plt.subplot(4,1,4)
plt.imshow(rel4.T,cmap='seismic_r',aspect='auto',interpolation='none', vmin=-0.5, vmax=0.50)
plt.ylabel('residue in S4')
ax = plt.gca();
ax.set_yticks(np.arange(0, 23, 1));
ax.set_yticklabels(np.arange(198, 221, 1));

plt.xlabel('every 100th frame')
plt.xticks(np.arange(0, 50, 5))
colorbar = plt.colorbar(cax = fig.add_axes([0.95, 0.12, 0.03, 0.76])) # left bot wid height
plt.savefig(key+'_byChainRel.eps',bbox_inches='tight')
