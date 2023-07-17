import numpy as np
from nefertiti.functions.superimpose import superimpose_array
import sys
from tqdm import tqdm

MAX_MATRIX = 20000
#MAX_MATRIX = 4000 ###

arr = np.load(sys.argv[1]) # $m-aa-fit-clust0.2.npy
#arr = arr[::5] ###
threshold = float(sys.argv[2])
big_cluster_file = sys.argv[3]
prmsd_file = sys.argv[4]

unclustered_mask = np.ones(len(arr), bool)

np.random.seed(0)
remain = arr

clusters = []
while len(remain) > MAX_MATRIX or len(clusters) == 0:
    print(f"Remaining structures: {len(remain)}")
    contest = remain
    indices = np.arange((len(remain)))
    for n in tqdm(range(10)):
        counts = np.zeros(len(contest), int)
        sample = remain[np.random.randint(len(remain), size=200)]
        sample2 = sample
        if len(contest) > 1000:
            sample2 = tqdm(sample)
        for s in sample2:
            _, rmsd = superimpose_array(contest, s)
            counts += (rmsd <= threshold)
        keep = np.argsort(-counts)[:int(len(contest)/2+0.5)]
        contest = contest[keep]
        indices = indices[keep]
    best = indices[0]
    _, rmsd = superimpose_array(remain, remain[best])
    in_cluster = (rmsd <= threshold)
    unclustered = np.where(unclustered_mask)[0]
    print(f"Structure {unclustered[best]+1}, cluster size {in_cluster.sum()}")
    clusters.append((unclustered[best], unclustered[in_cluster]))
    unclustered_mask[unclustered_mask] = ~in_cluster
    remain = arr[unclustered_mask]

with open(big_cluster_file, "w") as f:
    for cnr, (heart, cluster) in enumerate(clusters):
        print(f"Cluster {cnr+1} -> {heart+1}", end=" ", file=f)
        for ind in cluster:
            print(ind+1, end=" ", file=f)
        print("", file=f)    
    
mat = np.zeros((len(remain), len(remain)), np.float32)
for snr in tqdm(range(len(remain)-1)):
    _, rmsd = superimpose_array(remain[snr+1:], remain[snr])    
    mat[snr, snr+1:] = rmsd
    mat[snr+1:, snr] = rmsd

np.save(prmsd_file, mat)