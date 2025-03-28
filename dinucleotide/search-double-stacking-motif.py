import glob
import sys
import numpy as np
from numpy.linalg import norm, svd
from scipy.spatial.distance import pdist, squareform
from nefertiti.functions.superimpose import superimpose_array

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


find_double_stacking_out = sys.argv[1]
lineno = int(sys.argv[2])  # line in find_double_stacking.out
find_double_stacking = [l.split() for l in open(find_double_stacking_out).readlines()]
motif = find_double_stacking[lineno - 1]
print(motif)

refe_pdb = np.load(f"../database/pdbnpy/{motif[0]}.npy")
refe_pdb = refe_pdb[refe_pdb["model"] == 1]
refe_nuc = refe_pdb[refe_pdb["chain"] == motif[1].encode()]
refe_prot = refe_pdb[refe_pdb["chain"] == motif[6].encode()]

desc = calc_stacking_descriptors(motif[7], refe_prot)
desc.update(calc_stacking_descriptors(motif[9], refe_prot))
refe_desc1 = desc[int(motif[8])]
refe_desc2 = desc[int(motif[10])]
refe_desc = np.concatenate((refe_desc1, refe_desc2))
refe_desc_no_ca = np.concatenate((refe_desc1[1:], refe_desc2[1:]))

cbdis = norm(refe_desc[5] - refe_desc[1])
cgdis = norm(refe_desc[6] - refe_desc[2])

pdbs = glob.glob("../database/pdbnpy/*.npy")  # 2443 PDBs
###pdbs = ["3sxl_dom1.npy", "3sxl_dom2.npy"]  ###
###pdbs = ["/tmp/6y6e.npy"] ###
for pdbf in sorted(pdbs):
    #  code = pdbf  ###
    code = pdbf[pdbf.rindex("/") + 1 : pdbf.rindex(".")]
    # if code not in ("1B7F", "3TRZ", "3TS0", "3TS2", "4QQB", "5UDZ", "5YTT", "6A6J"):
    #    continue
    pdb = np.load(pdbf)
    pdb = pdb[pdb["model"] == 1]
    chains = np.unique(pdb["chain"])
    for chain in chains:
        pdbc = pdb[pdb["chain"] == chain]
        desc = {}
        for resname in ("PHE", "TYR", "TRP", "HIS", "ARG"):
            desc.update(calc_stacking_descriptors(resname, pdbc))
        if not desc:
            continue
        resids = list(desc.keys())
        resnames = {}
        for atom in pdbc:
            resid = atom["resid"]
            resnames[resid] = atom["resname"]
        alldesc = np.stack([desc[r] for r in resids], axis=0)

        coor_cb = alldesc[:, 1]
        coor_cg = alldesc[:, 2]
        ccbdist = squareform(pdist(coor_cb))
        ccgdist = squareform(pdist(coor_cg))

        MARGIN = 1.5
        mask = ((ccbdist - cbdis) ** 2 < MARGIN) & ((ccgdist - cgdis) ** 2 < MARGIN)
        p1, p2 = np.where(mask)
        if not p1.nbytes:
            continue

        """
        # superimpose, without CA
        pairdesc = np.empty((len(p1), 6, 3))
        pairdesc[:, :3] = alldesc[p1, 1:]
        pairdesc[:, 3:] = alldesc[p2, 1:]
        rotmat, rmsd = superimpose_array(pairdesc, refe_desc_no_ca)
        """

        # superimpose, with CA
        pairdesc = np.empty((len(p1), 8, 3))
        pairdesc[:, :4] = alldesc[p1]
        pairdesc[:, 4:] = alldesc[p2]
        rotmat, rmsd = superimpose_array(pairdesc, refe_desc)

        for n in rmsd.argsort():
            if rmsd[n] > 0.2:
                break

            resid1, resid2 = resids[p1[n]], resids[p2[n]]
            print(
                code,
                chain.decode(),
                resnames[resid1].decode(),
                resid1,
                resnames[resid2].decode(),
                resid2,
                "%.3f" % rmsd[n],
            )
        # sys.exit(0)
    # break
