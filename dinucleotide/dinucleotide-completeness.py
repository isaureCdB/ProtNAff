"""For each dinucleotide XY (and YZ) in XYZ-aa-fit-clust0.2.npy, 
 find the closest dincleotide in [XY?/?XY]-aa-fit-clust0.2.npy (and [YZ?/?YZ]-aa-fit-clust0.2.npy)

"""

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
POS = int(sys.argv[2])
assert POS in (0, 1)
DIHEXMOD = int(sys.argv[3])
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
for motif in motifs:
    # print(motif)
    frag = frags[motif]
    f = f"database/trilib/{motif}-aa-fit-clust0.2"
    clusters = read_clusterfile(f)
    f = f"database/trilib/{motif}-aa-fit-clust0.2.npy"
    coordinates[motif] = np.load(f)
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

scoordinates = {}
for motif in motifs:
    for pos in 0, 1:
        smotif = motif[pos : pos + 2]
        natoms = 22 * smotif.count("A") + 20 * smotif.count("C")
        if pos == 0:
            subc = coordinates[motif][:, :natoms]
        else:
            subc = coordinates[motif][:, -natoms:]
        subc = subc - subc.mean(axis=1)[:, None, :]
        scoordinates[motif, pos] = subc

motif = MOTIF
nfrag = len(coordinates[motif])
nsingletons = len([s for s in pdb_codes[motif] if s is not None])
print(f"{motif}: {nsingletons}/{nfrag} singletons", file=sys.stderr)
pos = POS
smotif = motif[pos : pos + 2]
natoms = 22 * smotif.count("A") + 20 * smotif.count("C")
# print(motif, smotif, natoms)
curr_pdb_codes = pdb_codes[motif]
for ind, coor in enumerate(scoordinates[motif, pos]):
    if ind % 32 != DIHEXMOD - 1:
        continue
    pdb_code = curr_pdb_codes[ind]
    if pdb_code is None:
        print(motif, pos, ind, 0)
        continue

    residuals1 = (coor * coor).sum()

    lowest_rmsd = 99999
    lowest_rmsd_motif = None
    lowest_rmsd_pos = None
    lowest_rmsd_clust = None

    for other_motif in motifs:
        mask = np.array([c == pdb_code for c in pdb_codes[other_motif]], dtype=bool)
        for other_pos in 0, 1:
            other_smotif = other_motif[other_pos : other_pos + 2]
            if other_smotif != smotif:
                continue
            # print("  ", other_motif, other_pos)
            # print("   MASK", mask.sum())
            other_coors = scoordinates[other_motif, other_pos]

            residuals2 = np.einsum("ijk,ijk->i", other_coors, other_coors)
            covar = np.einsum("ijk,jl->ikl", other_coors, coor)

            sd = superimpose_from_covar(covar, residuals1, residuals2)
            sd[mask] = np.inf
            argmin = sd.argmin()

            curr_lowest_rmsd = np.sqrt((sd[argmin] / natoms))
            if curr_lowest_rmsd < lowest_rmsd:
                lowest_rmsd = curr_lowest_rmsd
                lowest_rmsd_motif = other_motif
                lowest_rmsd_pos = other_pos
                lowest_rmsd_clust = argmin
    print(
        motif,
        pos,
        ind,
        "{:.3f}".format(lowest_rmsd),
        lowest_rmsd_motif,
        lowest_rmsd_pos,
        lowest_rmsd_clust,
    )
