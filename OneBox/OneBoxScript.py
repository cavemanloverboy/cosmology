import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree as Tree


# Load Data
# rpos is real space
# pos is redshift space
data = np.load('OneBox.npz')
import pdb; pdb.set_trace();
rpos, zpos = data.values()
boxsize = 1e3 # Mpc/h


# Generate randoms
Nrand = 10**6
rands = np.random.uniform(size=(Nrand,3))*boxsize

# Build trees
rtree = Tree(rpos, leafsize=4, compact_nodes=True, balanced_tree=True,boxsize=boxsize)
ztree = Tree(zpos, leafsize=4, compact_nodes=True, balanced_tree=True,boxsize=boxsize)


# Query Tree
# Gathers random-to-data distances (and ids)
k = [1]
r, ids   = rtree.query(rands, k=k)
zr, zids = ztree.query(rands, k=k)

# print Means
print(r.mean())
print(zr.mean())

# Compute pCDF and sort distances
cdf = np.arange(1, len(zr)+1)/len(zr)
pcdf = np.minimum(cdf, 1-cdf)
r  = np.sort( r[:,0]) # Quicksort is default algorithm , n log n
zr = np.sort(zr[:,0])


plt.figure()
plt.loglog( r, pcdf, label='real')
plt.loglog(zr, pcdf, label='redshift')
plt.legend(title='space')
plt.grid()
plt.xlabel('Distance (Mpc/h)')
plt.ylabel('Peaked CDF')
plt.ylim(1e-3,1)
plt.xlim(1,30)
plt.show()
