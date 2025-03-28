import json
import numpy as np

from numpy.linalg import svd, det


def superimpose_from_covar(covar, residuals1, residuals2):
    v, s, wt = svd(covar)
    reflect = det(v) * det(wt)
    s[:, -1] *= reflect
    sd = (residuals1 + residuals2) - 2 * s.sum(axis=1)
    sd = np.maximum(sd, 0)
    return sd


# print("Load fragment.json")
with open("database/fragments.json") as f:
    frags = json.load(f)

motifs = []
for n1 in ("A", "C"):
    for n2 in ("A", "C"):
        for n3 in ("A", "C"):
            motifs.append(n1 + n2 + n3)


def read_clusterfile(clusterfile):
    clusters = []
    for lnr, l in enumerate(open(clusterfile)):
        ll = l.split()
        assert ll[0].lower() == "cluster"
        assert ll[1] == str(lnr + 1)
        assert ll[2] == "->" or ll[2] == "=>"
        c = [int(lll) for lll in ll[3:]]
        clusters.append(c)
    return clusters


tot_nsingletons0 = 0
tot_nsingletons = 0
tot_clus = 0
for motif in motifs:
    # print(motif)
    frag = frags[motif]
    f = f"database/trilib/{motif}-aa-fit-clust0.2"
    clusters = read_clusterfile(f)
    f = f"database/trilib/{motif}-aa-fit-clust0.2.npy"
    curr_pdb_codes = []
    for clusnr, clus in enumerate(clusters):
        tot_clus += 1
        curr_pdb_code_set = set()
        for fr in clus:
            d = frag[str(fr)]
            struc = d["structure"][:4]
            curr_pdb_code_set.add(struc)
        if len(curr_pdb_code_set) == 1:
            tot_nsingletons += 1
        if len(clus) == 1:
            tot_nsingletons0 += 1

print(tot_clus)
print(tot_nsingletons)
print(tot_nsingletons0)
