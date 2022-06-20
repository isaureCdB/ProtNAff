# Functions taken from the filtering-clustering notebook

import numpy as np
from rmsdlib import multifit
from scipy.spatial.distance import pdist, squareform

def reassign(structures, clustids, chunksize):
    """Reassigns structures to the closest cluster center.
    Cluster centers are supplied in the form of cluster IDs, i.e.
    a list of integer IDs where each ID is the index of the cluster center.
    For example, [0, 3, 5] defines three cluster centers, the 1st 4th and 5th structure

    chunksize: number of structures to put in a chunk
        This is an implementation detail that only affects the speed, not the result
    """
    if len(structures.shape) == 3:
        assert structures.shape[2] == 3
        structures = structures.reshape(structures.shape[0], structures.shape[1]*3)
    if len(structures.shape) == 2:
        assert structures.shape[1] % 3 == 0

    clusters = {a:[a] for a in clustids}
    clus = structures[clustids]
    for n in range(0, len(structures), chunksize):
        chunk = structures[n:n+chunksize]
        d = chunk[:, np.newaxis, :] - clus[np.newaxis, :, :]
        inter_sd = np.einsum("...ij,...ij->...i", d,d)
        best = np.argmin(inter_sd, axis=1)
        for nn in range(len(chunk)):
            bestclust = clustids[best[nn]]
            if bestclust == (n+nn):
                continue
            clusters[bestclust].append(n+nn)
    return clusters

def cluster(structures, threshold, **kwargs):
    """Clusters structures using an RMSD threshold
        First, calculates the full pairwise RMSD matrix.
        The maximum number of structures is fixed at 20 000 (400 million pairwise RMSDs)
        
        The structure with the largest number of other structures within the RMSD
         threshold becomes a cluster.
        That structure and all of those other structures are removed.
        The process is repeated until no structures are left.

    This algorithm is used in ATTRACT and in the ProtNAff pipeline 
    (./create_frag_library/cluster-npy.py)
    """
    if len(structures.shape) == 3:
        assert structures.shape[2] == 3
        structures = structures.reshape(structures.shape[0], structures.shape[1]*3)
    if len(structures.shape) == 2:
        assert structures.shape[1] % 3 == 0

    if len(structures) > 20000:
        raise ValueError("Number of structures is too large, maximum is 20 000")

    dist = pdist(structures, 'sqeuclidean')
    d = squareform(dist)
    lim = threshold*threshold*structures.shape[1]/3
    d2 = d<lim
    del d

    clustids = []
    clustered = 0
    while clustered < len(structures):
        neigh = d2.sum(axis=0)
        cluster_center = neigh.argmax()
        clustids.append(cluster_center)
        cluster = np.where(d2[cluster_center])[0]
        for cs in cluster:
            d2[cs,:] = False
            d2[:, cs] = False
        clustered += len(cluster)
    return clustids

def fastcluster(structures, threshold, chunksize, existing=[], **kwargs):
    """Clusters structures using an RMSD threshold
        First structure becomes a cluster,
         second structure only if it doesn't cluster with the first, etc.

        structures: 2D numpy array, second dimension = 3 * natoms
          structures must already have been fitted!
        threshold: RMSD threshold (A)
        existing: a list of existing cluster IDs
        chunksize: number of structures to put in a chunk
          This is an implementation detail that only affects the speed, not the result
    
    This algorithm is used in ATTRACT and in the ProtNAff pipeline 
    (./create_frag_library/fastcluster_npy.py)
    """
    if len(structures.shape) == 3:
        assert structures.shape[2] == 3
        structures = structures.reshape(structures.shape[0], structures.shape[1]*3)
    if len(structures.shape) == 2:
        assert structures.shape[1] % 3 == 0
    natoms = structures.shape[1]/3
    # threshold2 = sum-of-sd threshold = (RMSD threshold **2) * natoms
    threshold2 = threshold**2 * natoms

    nclus = 1
    if len(existing):
        clus_space = max(len(existing), 100)
        clus = np.zeros((clus_space, structures.shape[1]))
        clus[:len(existing)] = structures[existing]
        clustids = []
        clustids[:] = existing[:]
    else:
        clus_space = 100
        clus = np.zeros((clus_space, structures.shape[1]))
        clus[:1] = structures[:1]
        clustids = [0]
    for n in range(1, len(structures), chunksize):
        #print("{0}/{1}".format(n, len(structures), file=sys.stderr)
        #sys.stderr.flush()
        chunk = structures[n:n+chunksize]
        d = chunk[:, np.newaxis, :] - clus[np.newaxis, :, :]
        inter_sd = np.einsum("...ij,...ij->...i", d, d)
        #close_inter is a 2D Boolean matrix:
        #  True  (1): chunk[i] is close to (within RMSD threshold of) clus[j]
        #  False (0): chunk[i] is not close to clus[j]
        close_inter = (inter_sd < threshold2)
        # newclustered contains all structures in the chunk that *don't* cluster with an existing cluster
        newclustered = []
        for chunk_index, closest_inter in enumerate(np.argmax(close_inter,axis=1)):
            # closest_inter contains the *first* index of close_inter
            #   with the highest value of close_inter
            # We are interested in the case where close_inter is all False (=> new cluster)
            # In that case, the highest value of close_inter is False, and closest_inter is 0
            # If close_inter is *not* all False (=> existing cluster), one of these conditions is False
            if closest_inter == 0 and close_inter[chunk_index, 0] == False:
                newclustered.append(chunk_index)

        if len(newclustered):
            # Now we have newclustered: the *chunk* index of all structures in the chunk that will be in new clusters
            # Now we want to cluster them among themselves, and add the *structure* id of each new cluster
            chunk_newclustered = chunk[newclustered]
            d = chunk_newclustered[:, np.newaxis, :] - chunk_newclustered[np.newaxis, :, :]
            intra_sd = np.einsum("...ij,...ij->...i", d, d)
            close_intra = (intra_sd < threshold2)

            # set all upper-triangular indices to False
            close_intra[np.triu_indices(len(chunk_newclustered))] = 0
            for nn in range(len(chunk_newclustered)):
                # same logic as for closest_inter;
                #  except that we don't have the chunk index, but the chunk_newclustered index (nn)
                #  and, since we modify close_intra in the "else" clause, argmax is computed later
                closest_intra = np.argmax(close_intra[nn])
                if closest_intra == 0 and close_intra[nn, 0] == False:
                    chunk_index = newclustered[nn]
                    # if clus is full, re-allocate it as a 50 % larger array
                    if nclus == clus_space:
                        clus_space = int(clus_space*1.5)
                        clus_old = clus
                        clus = np.zeros((clus_space, structures.shape[1]))
                        clus[:nclus] = clus_old
                    clus[nclus] = chunk[chunk_index]
                    clustids.append(n+chunk_index)
                    nclus += 1
                else:  # in addition, if we aren't a new cluster, being close to us doesn't matter
                    close_intra[:, nn] = False

    return clustids

def hierarchical_cluster(structures,threshold, secondary_threshold, chunksize, **kwargs):
    """Clusters structures using an RMSD threshold
    Does an initial clustering using fastcluster() with a secondary threshold
     then the resulting clusters are clustered using cluster()
    
    This is very close to the approach of the ProtNAff pipeline
    (./create_frag_library/clusterfrag_npy.sh)
    """
    if secondary_threshold <= threshold:
        raise ValueError("Secondary threshold must be higher than clustering threshold")
    print("Fast cluster", file=sys.stderr)
    initial_clustids = fastcluster(structures,secondary_threshold,chunksize=chunksize)
    print("Reassign", file=sys.stderr)
    initial_clusters = reassign(structures, initial_clustids, chunksize=chunksize)
    all_rough_clustids = []
    pos = 0
    for initial_cluster in initial_clusters.values():
        pos += 1
        print("{}/{}, size {}".format(pos, len(initial_clusters), len(initial_cluster)), file=sys.stderr)
        if len(initial_cluster) > 20000:
            raise ValueError("Initial cluster is too large, maximum is 20 000. Decrease secondary RMSD threshold.")
        initial_cluster = np.array(initial_cluster)
        initial_cluster_struc = structures[initial_cluster]
        subclustids0 = cluster(initial_cluster_struc, threshold)
        subclustids = initial_cluster[subclustids0]
        all_rough_clustids.append(subclustids)
    rough_clustids = np.concatenate(all_rough_clustids)
    rough_clust_struc = structures[rough_clustids]
    print("Recluster 1", file=sys.stderr)
    if len(rough_clustids) > 20000:
        clustids0 = fastcluster(rough_clust_struc, threshold, chunksize=chunksize)
    else:
        clustids0 = cluster(rough_clust_struc, threshold)
    fine_clustids = rough_clustids[clustids0].tolist()
    print("Recluster 2", file=sys.stderr)
    final_clustids = fastcluster(structures, threshold, chunksize=chunksize, existing=fine_clustids)
    return final_clustids    

def fit_multi_npy(a, ref):
    rotation, translation, RMSD = multifit(a, ref)
    rot = np.transpose(rotation, axes=(0,2,1))
    COM = a.sum(axis=1)/a.shape[1]
    centered = a - COM[:,None,:]
    rotated = np.einsum('...ij,...jk->...ik',centered,rot)
    fitted = rotated + COM[:,None,:]
    translated = fitted - translation[:,None,:]
    return translated, RMSD

if __name__ == "__main__":
    import sys
    struc = np.load(sys.argv[1])
    assert struc.ndim == 3, struc.shape
    assert struc.shape[2] == 3, struc.shape
    struc_fit_first = fit_multi_npy(struc, struc[0])[0]
    threshold = float(sys.argv[2])
    secondary_threshold = float(sys.argv[3])
    chunksize = 200
    final_clustids = hierarchical_cluster(struc_fit_first, threshold, secondary_threshold, chunksize)
    final_clusters = reassign(struc_fit_first, final_clustids, chunksize)

    for cnr, clustid in enumerate(final_clustids):
        clust = final_clusters[clustid]
        assert clustid == clust[0]
        print("Cluster {} -> {}".format(cnr+1, " ".join([str(n+1) for n in clust])))
    
    '''
    # verify
    print("Verify", file=sys.stderr)
    clusts = struc[final_clustids]
    for snr, s in enumerate(struc): 
        print("{}/{}".format(snr+1, len(struc)), file=sys.stderr)
        rmsds = fit_multi_npy(clusts, s)[1]
        assert rmsds.min() <= threshold, (snr, rmsds.min())
    '''

    '''
    # Verify hierarchique
    h = "41 52 459 313 572 893 927 565 469 552 554 598 892"
    h_clustids = np.array([int(ss) for ss in h.split()]) - 1
    print(h_clustids)
    clusts = struc[h_clustids]
    for snr, s in enumerate(struc): 
        rmsds = fit_multi_npy(clusts, s)[1]
        assert rmsds.min() <= threshold, (snr, rmsds.min())
    '''