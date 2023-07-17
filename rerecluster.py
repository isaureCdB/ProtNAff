import numpy as np
import sys
from nefertiti.functions.superimpose import superimpose_array
from tqdm import tqdm

arr = np.load(sys.argv[1])  # $m-aa-fit-clust0.2.npy
big_cluster_list = sys.argv[2]
prmsd_file = sys.argv[3]
original_rmsd = float(sys.argv[4])
new_rmsd = float(sys.argv[5])
output_cluster_list_all = sys.argv[6]
output_cluster_list_closest = sys.argv[7]
output_cluster_closest_rmsd = sys.argv[8]

clusters = []
for lnr, l in enumerate(open(big_cluster_list)):
    ll = l.split()
    assert ll[0] == "Cluster"
    assert ll[1] == str(lnr+1)
    assert ll[2] == "->"
    heart = int(ll[3])
    c = np.array(list(set([int(lll) for lll in ll[3:]])))
    clusters.append((heart, c))

prmsd = np.load(prmsd_file)
prmsd_mask = (prmsd <= new_rmsd)

unclustered_mask = np.ones(len(arr), bool)
for heart, c in clusters:
    unclustered_mask[c-1] = 0

remain = arr[unclustered_mask]
assert len(remain) == len(prmsd) == prmsd.shape[1]
remain_indices = np.where(unclustered_mask)[0]

print("Cluster pre-calculated pairwise RMSD matrix")
remain_clusters = []  # index in prmsd, starting from 0
remain_clustered = 0
progress = tqdm(total=len(remain))
while remain_clustered < len(remain):
    nb = prmsd_mask.sum(axis=0)
    best = nb.argmax()
    c = np.where(prmsd_mask[best])[0]
    assert len(c)
    assert len(c) == nb[best]
    '''
    # check
    ch = remain[best]
    _, rmsds = superimpose_array(remain[c], ch)
    assert rmsds.max() < new_rmsd
    # /check
    '''
    if len(c) == 1:
        diag = np.diag(prmsd_mask)
        singletons = np.where(diag)[0]
        for vnr in range(len(singletons)):
            remain_clusters.append(singletons[vnr])
        remain_clustered += len(singletons)
        progress.update(len(singletons))
        break
    remain_clusters.append(best)
    remain_clustered += len(c)
    prmsd_mask[c, :] = 0
    prmsd_mask[:, c] = 0
    progress.update(len(c))
progress.close()

reclusters = remain_indices[remain_clusters].tolist()

print("Re-cluster previous big clusters")
t = 2 * original_rmsd + new_rmsd
for heart, c in tqdm(clusters):
    existing_cluster_struc = arr[reclusters]
    _, heart_rmsd = superimpose_array(existing_cluster_struc, arr[heart]) 
    big_cluster_struc = arr[c-1]
    keep = np.ones(len(big_cluster_struc), bool)
    # Find the big cluster members that are close to an existing cluster
    for cnr, c2 in enumerate(tqdm(existing_cluster_struc[heart_rmsd < t])):
        _, rmsds = superimpose_array(big_cluster_struc, c2)
        low = (rmsds<new_rmsd)
        to_remove = keep & low
        keep[to_remove] = 0
    inds = (c-1)[keep]
    keep_struc = big_cluster_struc[keep]
    mat = np.ones((len(keep_struc), len(keep_struc)), bool)
    for snr in tqdm(range(len(keep_struc)-1)):
        _, rmsd = superimpose_array(keep_struc[snr+1:], keep_struc[snr])
        low = (rmsd < new_rmsd)
        mat[snr, snr+1:] = low
        mat[snr+1:, snr] = low
    
    reclustered = 0
    progress = tqdm(total=len(mat))
    while reclustered < len(mat):
        nb = mat.sum(axis=0)
        best = nb.argmax()
        c = np.where(mat[best])[0]
        assert len(c)
        assert len(c) == nb[best]
        if len(c) == 1:
            diag = np.diag(mat)
            singletons = np.where(diag)[0]
            for vnr in range(len(singletons)):
                best = singletons[vnr]
                assert np.all(arr[inds[best]] == keep_struc[best])
                reclusters.append(inds[best])
            reclustered += len(singletons)
            progress.update(len(singletons))
            break
        assert np.all(arr[inds[best]] == keep_struc[best])
        reclusters.append(inds[best])
        reclustered += len(c)
        mat[c, :] = 0
        mat[:, c] = 0
        progress.update(len(c))
    progress.close()

print("Assign all structures to the closest cluster and to all clusters within the threshold")
cluster_struc = arr[reclusters]
hearts = set(reclusters)
clustering_all = {h: [h] for h in hearts}
clustering_closest = {h: [h] for h in hearts}
for struc_nr, struc in enumerate(tqdm(arr)):
    _, rmsds = superimpose_array(cluster_struc, struc)
    if struc_nr in hearts:
        pos = reclusters.index(struc_nr)
        assert rmsds[pos] < 0.0001
        continue
    best = rmsds.argmin()
    assert rmsds[best] < new_rmsd, struc_nr + 1
    clustering_closest[reclusters[best]].append(struc_nr)
    for good in np.where(rmsds < new_rmsd)[0]:
        clustering_all[reclusters[good]].append(struc_nr)

def write_cluster(clusterfile, clusters):
    with open(clusterfile, "w") as f:
        for clust_nr, cluster in enumerate(clusters):
            print(f"Cluster {clust_nr+1} -> ", end=" ", file=f)
            for ind in cluster:
                print(ind+1, end=" ", file=f)
            print("", file=f)    

clusters = sorted(clustering_all.items(), key = lambda item: -999999 * len(item[1]) + item[0])
clustering_all = [clustering_all[k] for k,v in clusters]
clustering_closest = [clustering_closest[k] for k,v in clusters]

print("Write out clusters")
write_cluster(output_cluster_list_all, clustering_all)
write_cluster(output_cluster_list_closest, clustering_closest)

print("Analyze closest cluster for each cluster")
closest_rmsd = []
closest_cluster = []
cluster_struc = arr[[it[0] for it in clusters]]
for struc_nr, struc in enumerate(tqdm(cluster_struc)):
    _, rmsds = superimpose_array(cluster_struc, struc)
    rmsds[struc_nr] = 9999
    closest_rmsd.append(rmsds.min())
    closest_cluster.append(rmsds.argmin())
closest_rmsd = np.array(closest_rmsd)
closest_cluster = np.array(closest_cluster)

with open(output_cluster_closest_rmsd, "w") as f:
    for rmsd, clust in zip(closest_rmsd, closest_cluster):
        print("{:.3f} {:d}".format(rmsd, clust+1), file=f)
