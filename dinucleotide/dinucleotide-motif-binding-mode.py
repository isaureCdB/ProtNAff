"""For a dinucleotide, defined by the first two nucleotides in a 0.2A trinucleotide cluster,
 find all dinucleotides in any 0.2A trinucleotide cluster within 1A.
 Once identified, apply the rotation/translation matrices to the entire RNA and to the protein interface.
"""

import numpy as np

import sys

from nefertiti.functions.superimpose import superimpose_array, superimpose
from nefertiti.functions.parse_pdb import parse_pdb
from nefertiti.functions.write_pdb import write_pdb

def mutseq(seq):
    return seq.decode().replace("U", "C").replace("G", "A")


refe_motif = sys.argv[1]
refe_cluster = int(sys.argv[2])
output_dir = sys.argv[3]
pdb_list = None
if len(sys.argv) > 4:
    pdb_list = []
    pdb_list_file = sys.argv[4]
    for l in open(pdb_list_file):
        l = l.strip()
        if len(l):
            assert len(l) == 4, l
            pdb_list.append(l)

all_frags = np.load("database/fragments_clust.npy")
all_seqs = [mutseq(frag["seq"]) for frag in all_frags]
frags = {}

motifs = []
for n1 in ("A", "C"):
    for n2 in ("A", "C"):
        for n3 in ("A", "C"):
            motifs.append(n1 + n2 + n3)

coordinates = {}
frags = {}
tot_frags = {}
for motif in motifs:
    f = f"database/trilib/{motif}-aa-fit-clust0.2.npy"
    coordinates[motif] = np.load(f)
    frags0 = all_frags[[s == motif for s in all_seqs]]
    tot_frags[motif] = frags0
    frags1 = frags0[frags0["clust0.2A_center"] == 1]
    frags2 = frags1[frags1["clust0.2A"].argsort()]
    assert np.all(np.equal(frags2["clust0.2A"], np.arange(len(frags2)).astype(int) + 1))
    frags[motif] = frags2
    assert len(coordinates[motif]) == len(frags[motif]), (
        motif,
        len(coordinates[motif]),
        len(frags[motif]),
    )

scoordinates = {}
scoormean = {}
for motif in motifs:
    smotif = motif[:2]
    natoms = 22 * smotif.count("A") + 20 * smotif.count("C")
    subc = coordinates[motif][:, :natoms]
    scoormean[motif] = subc.mean(axis=1)[:, None, :]
    scoordinates[motif] = subc - scoormean[motif]

natoms = 22 * smotif.count("A") + 20 * smotif.count("C")

refe_coor = scoordinates[refe_motif][refe_cluster - 1]
refe_mean = scoormean[refe_motif][refe_cluster - 1]

close_frags = []
for other_motif in motifs:
    if other_motif[:2] != refe_motif[:2]:
        continue
    other_coors = scoordinates[other_motif]

    rotmat, rmsd = superimpose_array(other_coors, refe_coor)
    mask = rmsd <= 1
    if mask.sum():
        close_clust = frags[other_motif][mask]["clust0.2A"]
        mask2 = np.isin(tot_frags[other_motif]["clust0.2A"], close_clust)
        close_frags.append(tot_frags[other_motif][mask2])

close_frags = np.concatenate(close_frags)

for n, close_frag in enumerate(close_frags):
    if close_frag["model"] != 1:
        continue
    struc, chain = close_frag["structure"].decode(), close_frag["chain"].decode()
    if pdb_list and struc not in pdb_list:
        continue
    fname = f"database/cleanPDB/{struc}{chain}-1-iniparse-aa.pdb"
    try:
        with open(fname) as fp:
            pdbtxt = fp.read()
    except FileNotFoundError:
        print("skip ", struc, chain)
        continue
    pdb = parse_pdb(pdbtxt)
    pdb_coors = np.stack((pdb["x"], pdb["y"], pdb["z"]),axis=1)
    fit_mask = (pdb["resid"] == close_frag["indices"][0]) | (pdb["resid"] == close_frag["indices"][1])
    fit_mask = fit_mask & ~((pdb["resname"] == b'RG') & (pdb["name"] == b'N2'))
    assert fit_mask.sum() == len(refe_coor), (fit_mask.sum(), len(refe_coor))
    fit_coors = pdb_coors[fit_mask]
    com = fit_coors.mean(axis=0)
    rotmat, rmsd = superimpose(fit_coors - com, refe_coor)

    pdb_rcoors = (pdb_coors - com).dot(rotmat)
    pdb["x"] = pdb_rcoors[:, 0]
    pdb["y"] = pdb_rcoors[:, 1]
    pdb["z"] = pdb_rcoors[:, 2]


    ind = close_frag["indices"][0]

    print("write", struc, chain, ind)

    pdbtxt = write_pdb(pdb)
    fname = f"{output_dir}/{struc}_{chain}-{ind}-RNA.pdb"
    with open(fname, "w") as fp:
        fp.write(pdbtxt)

    pdbtxt = write_pdb(pdb[fit_mask])
    fname = f"{output_dir}/{struc}_{chain}-{ind}-frag.pdb"
    with open(fname, "w") as fp:
        fp.write(pdbtxt)

    fname = f"database/interface/{struc}-1.pdb"
    with open(fname) as fp:
        pdbtxt = fp.read()
    pdb = parse_pdb(pdbtxt)
    pdb_coors = np.stack((pdb["x"], pdb["y"], pdb["z"]),axis=1)
    pdb_rcoors = (pdb_coors - com).dot(rotmat)
    pdb["x"] = pdb_rcoors[:, 0]
    pdb["y"] = pdb_rcoors[:, 1]
    pdb["z"] = pdb_rcoors[:, 2]

    pdbtxt = write_pdb(pdb)
    fname = f"{output_dir}/{struc}_{chain}-{ind}-prot.pdb"
    with open(fname, "w") as fp:
        fp.write(pdbtxt)
