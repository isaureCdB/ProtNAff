import glob
import os
import sys
import numpy as np
from numpy.linalg import norm, svd

rings = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "ARG": ["CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "TRP": ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
}

ringsb = {k.encode(): np.array([vv.encode() for vv in v]) for k, v in rings.items()}


def calc_plane(res, rescoor):
    resname = res[0]["resname"]
    mask = np.isin(res["name"], ringsb[resname])
    coor = rescoor[mask]
    cavec = np.zeros(3)
    for n, atom in enumerate(res):
        name = atom["name"]
        if name == b"CA":
            cavec += rescoor[n]
        elif name == b"CB":
            cavec -= rescoor[n]

    coor = coor - coor.mean(axis=0)
    covar = coor.T.dot(coor)
    v, s, wt = svd(covar)
    plane = wt[2] / norm(wt[2])
    if plane.dot(cavec) > 0:
        plane *= -1
    return plane


def calc_planes(coor, ca, cb):
    coor = coor - coor.mean(axis=1)[:, None]
    covar = np.einsum("ijk,ijl->ikl", coor, coor)  # coor[i].T.dot(coor[i])
    v, s, wt = svd(covar)
    planes_directionless = wt[:, 2, :] / norm(wt[:, 2, :], axis=-1)[..., None]
    cavec = ca - cb
    flip = np.sign(np.einsum("ij,ij->i", planes_directionless, -cavec))
    return planes_directionless * flip[:, None]


def calc_stacking_descriptors(resname, chain_pdb):
    atoms = [b"CA", b"CB"] + [r.encode() for r in rings[resname]]
    aa_mask = chain_pdb["resname"] == resname.encode()
    aa0 = chain_pdb[aa_mask]

    atom_masks0 = []
    for atom in atoms:
        atom_mask0 = aa0["name"] == atom
        atom_masks0.append(atom_mask0)
    common = set(aa0[atom_masks0[0]]["resid"])
    for atom_mask0 in atom_masks0:
        common = common.intersection(set(aa0[atom_mask0]["resid"]))

    if not len(common):
        return {}
    common_mask = np.isin(aa0["resid"], list(common))
    aa = aa0[common_mask]
    aa_resids = None
    aa_coor0 = []
    for atom in atoms:
        mask = aa["name"] == atom
        aaa = aa[mask]
        if len(aaa) != len(common):
            print("ERR", file=sys.argv[1])
            return {}
        if atom == b"CA":
            aa_resids = aaa["resid"]
        aaa_coor = np.stack((aaa["x"], aaa["y"], aaa["z"]), axis=1)
        aa_coor0.append(aaa_coor)
    aa_coor = np.empty((len(common), len(atoms), 3))
    for n, aa_c0 in enumerate(aa_coor0):
        aa_coor[:, n, :] = aa_c0
    desc = np.empty((len(common), 4, 3))
    desc[:, :3] = aa_coor[:, :3]  # copy CA, CB and CG
    desc_plane = calc_planes(aa_coor[:, 2:], aa_coor[:, 0], aa_coor[:, 1])
    desc[:, 3] = desc_plane * 1.4 + aa_coor[:, 2]

    result = {}
    for x, y in zip(aa_resids, desc):
        result[int(x)] = y
    return result


curr_code = None
miss = False

find_double_stacking_out = sys.argv[1]
sorted_double_stacking_out = sys.argv[2]
output_descriptors_file = sys.argv[3]

find_double_stacking = [l.split() for l in open(find_double_stacking_out).readlines()]
find_double_stacking.sort()

refe_descs = []

outf = open(sorted_double_stacking_out, "w")

for motif in find_double_stacking:

    if motif[0] != curr_code:
        curr_code = motif[0]
        pdbf = f"../database/pdbnpy/{motif[0]}.npy"
        if not os.path.exists(pdbf):
            print("MISS", curr_code)
            miss = True
            continue
        miss = False
        refe_pdb = np.load(pdbf)
        refe_pdb = refe_pdb[refe_pdb["model"] == 1]
    elif miss:
        continue
    # refe_nuc = refe_pdb[refe_pdb["chain"] == motif[1].encode()]
    refe_prot = refe_pdb[refe_pdb["chain"] == motif[6].encode()]

    desc = calc_stacking_descriptors(motif[7], refe_prot)
    desc.update(calc_stacking_descriptors(motif[9], refe_prot))
    try:
        refe_desc1 = desc[int(motif[8])]
        refe_desc2 = desc[int(motif[10])]
    except KeyError:
        print("ERR", curr_code)
        continue
    print(" ".join(motif))
    print(" ".join(motif), file=outf)
    refe_descs.append((refe_desc1, refe_desc2))

print(len(refe_descs))
refe_descs1 = np.stack([d[0] for d in refe_descs], axis=0)
refe_descs2 = np.stack([d[1] for d in refe_descs], axis=0)
refe_descs = np.concatenate((refe_descs1, refe_descs2), axis=1)
np.save(output_descriptors_file, refe_descs)
