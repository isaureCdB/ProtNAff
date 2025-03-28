import sys
import numpy as np
from scipy.spatial.distance import squareform
from nefertiti.functions.superimpose import superimpose_array
from scipy.cluster.hierarchy import linkage, dendrogram

descriptors_file = sys.argv[1]
descriptors = np.load(descriptors_file)
descriptors_rmsd_outfile = sys.argv[2]

rmsds = []
for n, desc in enumerate(descriptors):
    print(n + 1, len(descriptors))
    desc_array = descriptors[n + 1 :]
    if not len(desc_array):
        break
    _, rmsd = superimpose_array(desc_array, desc)
    rmsds.append(rmsd)
rmsds_condensed = np.concatenate(rmsds)
rmsds_sq = squareform(rmsds_condensed)
np.save(descriptors_rmsd_outfile, rmsds_sq)
