# Script version of clustering.ipynb and colab/clustering.ipynb
# For quick testing/debugging

#####################################
# Fragment library definition (before loading)
# These are global variables
#####################################

motifs = ["AAA", "AAC", "ACA", "ACC", "CAA", "CAC", "CCA", "CCC"]
structures = None
fragments = None
coordinates = None
rev_fragindex = None  # to map a key in fragments dict to an index in the fragment array

#####################################
# Masks for filtering
#####################################

def mask_pure_ss(motif_fragments,structures, motif_rev_fragindex):
    """Returns a mask that indicates for each fragment 
if it is purely single-stranded or not"""
    ss_set = set(["L", "T", "S", "J", "B", "I"])
    
    reverse_mapping = {}
    for strucname, structure in structures.items():
        struc_rev_mapping = {}
        for chain, mapping in structure["mapping"].items():
            struc_rev_mapping[chain] = {v:k for k,v in mapping.items()}
        reverse_mapping[strucname] = struc_rev_mapping
    
    result = np.zeros(len(motif_rev_fragindex),bool)

    for fragkey, frag in motif_fragments.items():
        if fragkey not in motif_rev_fragindex:
            # fragment contains clashes
            continue
        fragind = motif_rev_fragindex[fragkey]
        structure = frag["structure"]
        chain = frag["chain"]
        rmap = reverse_mapping[structure]["chain_" + chain]
        res = [rmap[resid] for resid in frag["resid"]]
        ss_states = structures[structure]["ss"]["chain_" + chain]
        frag_ss_states = [ss_states["res_"+ r][0] for r in res]
        is_pure_ss = all([ss_state in ss_set for ss_state in frag_ss_states])
        result[fragind] = is_pure_ss

    return result


def mask_nmr(motif_fragments,structures, motif_rev_fragindex):
    """Returns a mask that indicates for each fragment 
if it comes from NMR or not"""
    result = np.zeros(len(motif_rev_fragindex),bool)

    for fragkey, frag in motif_fragments.items():
        if fragkey not in motif_rev_fragindex:
            # fragment contains clashes
            continue
        fragind = motif_rev_fragindex[fragkey]
        structure = frag["structure"]
        struc = structures[structure]
        is_nmr = (struc["method"].find("nmr") > -1)
        result[fragind] = is_nmr
    
    return result

######################################
# Mask-based filters
#####################################

def filter_pure_ss():
    result = {}
    for motif in fragments:
        frags = fragments[motif]
        mask = mask_pure_ss(frags, structures, rev_fragindex[motif])
        result[motif] = coordinates[motif][mask]
    return result

def filter_pure_ss_nmr():
    result = {}
    for motif in fragments:
        frags = fragments[motif]
        args = frags, structures, rev_fragindex[motif]
        mask1 = mask_pure_ss(*args)
        mask2 = mask_nmr(*args)
        mask = mask1 & mask2
        result[motif] = coordinates[motif][mask]
    return result

def filter_pure_ss_not_nmr():
    result = {}
    for motif in fragments:
        frags = fragments[motif]
        args = frags, structures, rev_fragindex[motif]
        mask1 = mask_pure_ss(*args)
        mask2 = mask_nmr(*args)
        mask = mask1 & (~mask2)
        result[motif] = coordinates[motif][mask]
    return result

def filter_none():
    return coordinates

filters = {
    "None": filter_none,
    "Pure single-stranded": filter_pure_ss,
    "Pure single-stranded, only NMR structures": filter_pure_ss_nmr,
    "Pure single-stranded, no NMR structures": filter_pure_ss_not_nmr,
}

######################################
# Clustering methods
#####################################

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
    chunksize = 200
    initial_clustids = fastcluster(structures,secondary_threshold,chunksize=chunksize)
    initial_clusters = reassign(structures, initial_clustids, chunksize=chunksize)
    all_rough_clustids = []
    for initial_cluster in initial_clusters.values():
        if len(initial_cluster) > 20000:
            raise ValueError("Initial cluster is too large, maximum is 20 000. Increase secondary RMSD threshold.")
        initial_cluster = np.array(initial_cluster)
        initial_cluster_struc = structures[initial_cluster]
        subclustids0 = cluster(initial_cluster_struc, threshold)
        subclustids = initial_cluster[subclustids0]
        all_rough_clustids.append(subclustids)
    rough_clustids = np.concatenate(all_rough_clustids)
    rough_clust_struc = structures[rough_clustids]
    if len(rough_clustids) > 20000:
        clustids0 = fastcluster(rough_clust_struc, threshold, chunksize=chunksize)
    else:
        clustids0 = cluster(rough_clust_struc, threshold)
    fine_clustids = rough_clustids[clustids0].tolist()
    final_clustids = fastcluster(structures, threshold, chunksize=chunksize, existing=fine_clustids)
    return final_clustids    

clustering_methods = {
    "Fast": fastcluster,
    "Full matrix": cluster,
    "Hierarchical": hierarchical_cluster,
}

#####################################
# User parameters
#####################################

filtering_method = "Pure single-stranded, no NMR structures"
clustering_method = "Hierarchical"
clustering_threshold = 1.0
secondary_threshold = 2.0
chunksize = 200  # affects only the speed, not the result

#####################################
# Parameter validation
#####################################

if filtering_method not in filters:
    raise ValueError("Filtering method '{}' not in {}".format(filtering_method, filters.keys()))
if clustering_method not in clustering_methods:
    raise ValueError("Clustering method '{}' not in {}".format(clustering_method, clustering_methods.keys()))
clustering_threshold = float(clustering_threshold)
if clustering_threshold <= 0:
    raise ValueError("Clustering threshold must be positive")
secondary_threshold = float(secondary_threshold)
chunksize = int(chunksize)
if chunksize <= 0:
    raise ValueError("Chunk size must be positive")

#####################################
# Imports
#####################################

# In colab: install scipy with pip

import os
import tempfile
import shutil
import json
import numpy as np
from scipy.spatial.distance import pdist, squareform
from urllib.request import urlretrieve

#####################################
# Loading fragment library
#####################################

fraglib_files = [
    "structures.json",

    "fragments_clust.json",

    "aa-fit.tgz",  
    # containing:
    #   <motif>-aa-fit.npy
    #   <motif>-aa.noclash
    #   for 8 motifs (AAA to CCC)
    # These come from the trilib/ directory.
    # The .npy contains only the non-clashing fragments

]

with tempfile.TemporaryDirectory() as TEMPDIR:
    for fraglib_file in fraglib_files:
        # dummy code: copy from ./example/full_test/for_colab to TEMPDIR
        # In the notebook, use urlretrieve
        src = "../example/full_test/for_colab/" + fraglib_file
        dest = os.path.join(TEMPDIR, fraglib_file)
        shutil.copyfile(src, dest)
        # /dummy code
    os.system("cd {}; tar xzf aa-fit.tgz".format(TEMPDIR))

    structures = json.load(open(os.path.join(TEMPDIR, "structures.json")))
    fragments = json.load(open(os.path.join(TEMPDIR, "fragments_clust.json")))
    coordinates = {}
    rev_fragindex = {}
    for motif in motifs:
        coorfile = "{}-aa-fit.npy".format(motif)
        coordinates[motif] = np.load(os.path.join(TEMPDIR, coorfile))
        noclash_file = "{}-aa.noclash".format(motif)
        motif_rev_fragindex = {}
        with open(os.path.join(TEMPDIR, noclash_file)) as f:
            for lnr, l in enumerate(f.readlines()):
                frag_key = l.strip()
                motif_rev_fragindex[frag_key] = lnr
            assert len(motif_rev_fragindex) == len(coordinates[motif]), (len(motif_rev_fragindex), noclash_file, len(coordinates[motif]), coorfile)
            rev_fragindex[motif] = motif_rev_fragindex


#####################################
# Run filtering and clustering
#####################################

filt = filters[filtering_method]
clustering_func = clustering_methods[clustering_method]
filtered_coordinates = filt()
clusters = {}
for motif in motifs:
    motif_coors = filtered_coordinates[motif]
    motif_clustids = clustering_func(
        motif_coors,
        threshold=clustering_threshold,
        secondary_threshold=secondary_threshold,
        chunksize=chunksize
    )
    motif_clusters = reassign(
        motif_coors,
        motif_clustids,
        chunksize=chunksize
    )
    clusters[motif] = motif_clusters

for motif in motifs:
    print(motif, len(clusters[motif]))