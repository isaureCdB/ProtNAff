import json
import numpy as np

import sys

from numpy.linalg import svd, det


def superimpose_from_covar(covar, residuals1, residuals2):
    v, s, wt = svd(covar)
    reflect = det(v) * det(wt)
    s[:, -1] *= reflect
    sd = (residuals1 + residuals2) - 2 * s.sum(axis=1)
    sd = np.maximum(sd, 0)
    return sd


MOTIF = sys.argv[1]
DIHEXMOD = int(sys.argv[2])
assert DIHEXMOD >= 1 and DIHEXMOD <= 32

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


pdb_codes = {}  # singleton PDB code, or None if not a singleton
coordinates = {}
all_residuals = {}
for motif in motifs:
    # print(motif)
    frag = frags[motif]
    f = f"database/trilib/{motif}-aa-fit-clust0.2"
    clusters = read_clusterfile(f)
    f = f"database/trilib/{motif}-aa-fit-clust0.2.npy"
    coordinates[motif] = np.load(f)
    all_residuals[motif] = np.einsum(
        "ijk,ijk->i", coordinates[motif], coordinates[motif]
    )
    assert len(coordinates[motif]) == len(clusters)
    curr_pdb_codes = []
    for clusnr, clus in enumerate(clusters):
        curr_pdb_code_set = set()
        for fr in clus:
            d = frag[str(fr)]
            struc = d["structure"][:4]
            curr_pdb_code_set.add(struc)
        pdb_code = curr_pdb_code_set.pop() if len(curr_pdb_code_set) == 1 else None
        curr_pdb_codes.append(pdb_code)
    pdb_codes[motif] = curr_pdb_codes

motif = MOTIF
nfrag = len(coordinates[motif])
nsingletons = len([s for s in pdb_codes[motif] if s is not None])
print(f"{motif}: {nsingletons}/{nfrag} singletons", file=sys.stderr)
natoms = coordinates[motif].shape[1]
curr_pdb_codes = pdb_codes[motif]

other_coors = coordinates[motif]
residuals2 = all_residuals[motif]
for ind, coor in enumerate(coordinates[motif]):
    if ind % 32 != DIHEXMOD - 1:
        continue
    pdb_code = curr_pdb_codes[ind]
    if pdb_code is None:
        print(motif, ind, 0)
        continue

    residuals1 = (coor * coor).sum()

    lowest_rmsd = 99999
    lowest_rmsd_clust = None

    mask = np.array([c == pdb_code for c in pdb_codes[motif]], dtype=bool)

    covar = np.einsum("ijk,jl->ikl", other_coors, coor)

    sd = superimpose_from_covar(covar, residuals1, residuals2)
    sd[mask] = np.inf
    argmin = sd.argmin()

    curr_lowest_rmsd = np.sqrt((sd[argmin] / natoms))
    if curr_lowest_rmsd < lowest_rmsd:
        lowest_rmsd = curr_lowest_rmsd
        lowest_rmsd_clust = argmin
    print(
        motif,
        ind,
        "{:.3f}".format(lowest_rmsd),
        lowest_rmsd_clust,
    )
