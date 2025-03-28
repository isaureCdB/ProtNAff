from bisect import bisect
import json
import sys
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram


descriptors_rmsd_file = sys.argv[1]
motifs_file = sys.argv[2]  # double-stacking-sorted.out
clustering_threshold = float(sys.argv[3])  # 0.5
out_dendrogram_png = sys.argv[4]
out_clustering_file = sys.argv[5]

motifs = [l.strip() for l in open(motifs_file)]
rmsds_sq = np.load(descriptors_rmsd_file)
assert rmsds_sq.shape == (len(motifs), len(motifs)), (rmsds_sq.shape, len(motifs))

rmsds_condensed = squareform(rmsds_sq)
linkage_matrix = linkage(rmsds_condensed, method="average")

order = dendrogram(
    linkage_matrix,
    truncate_mode=None,
    no_plot=True,
    show_contracted=True,
    distance_sort="descending",
)["leaves"]
ordered_motifs = [motifs[o] for o in order]

nclust = len(linkage_matrix) - bisect(linkage_matrix[:, 2], clustering_threshold)

d = dendrogram(
    linkage_matrix,
    truncate_mode="lastp",
    p=nclust,
    no_plot=True,
    show_leaf_counts=True,
    distance_sort="descending",
)
clustsizes = d["ivl"]

clustering = []
pos = 0
for label in clustsizes:
    if label[0] == "(":
        size = int(label[1:-1])
        clust = ordered_motifs[pos : pos + size]
        clustering.append(clust)
        pos += size
    else:
        o = int(label)
        assert o == order[pos]
        clustering.append([ordered_motifs[pos]])
        pos += 1
with open(out_clustering_file, "w") as fp:
    json.dump(clustering, fp, indent=2)

cluster_labels = []
for clustnr, clust in enumerate(clustering):
    codes = [c.split()[0] for c in clust]
    aa = [c.split()[7] + " " + c.split()[9] for c in clust]

    unique_codes, counts = np.unique(codes, return_counts=True)
    mode_code = unique_codes[counts.argmax()]
    _a, _b = np.unique(aa, return_counts=True)
    mode_aa = _a[_b.argmax()]
    xes = int(len(codes) / 5 + 0.5)
    cluster_labels.append(f"{clustnr+1:3d} {xes*'X'} {mode_code} {mode_aa}")

assert len(clustering) == len(d["leaves"])
cluster_label_dict = {k: v for k, v in zip(d["leaves"], cluster_labels)}

fig, ax = plt.subplots(figsize=(15, 15), dpi=500)
r = dendrogram(
    linkage_matrix,
    truncate_mode="lastp",
    orientation="left",
    leaf_font_size=2,
    p=nclust,
    distance_sort="descending",
    leaf_label_func=lambda i: cluster_label_dict[i],
)
ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=0.2))
ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(base=0.05))
ax.xaxis.set_tick_params(labelsize=2)
fig.savefig(out_dendrogram_png, bbox_inches="tight", pad_inches=0.1)
