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


_3_to_1 = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "MSE": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}


def get_seqsig(resid1, resid2, chain_ca):
    resids = chain_ca["resid"]
    s = []
    for resid in (resid1, resid2):
        ss = ["-"] * 15
        atoms = chain_ca[(resids >= resid - 7) & (resids <= resid + 7)]
        lastr = None
        for atom in atoms:
            if atom["resid"] != lastr:
                lastr = atom["resid"]
                pos = lastr - resid + 7
                # print(pos, lastr, resid)
                resname = atom["resname"].decode()
                ss[pos] = _3_to_1.get(resname, "X")
        s.append("".join(ss))
    return s[0] + "..." + s[1]


double_stacking_sorted_file = sys.argv[1]
lineno = int(sys.argv[2])  # line in find_double_stacking.out
motifs = [l.strip() for l in open(double_stacking_sorted_file).readlines()]
double_stacking_sorted_seqsig_file = sys.argv[3]
seqsigs = [l.strip() for l in open(double_stacking_sorted_seqsig_file).readlines()]
assert len(seqsigs) == len(motifs)

truepos_outfile = sys.argv[4]
falsepos_outfile = sys.argv[5]

motif_inds = {}
for ind, motif in enumerate(motifs):
    m = motif.split()
    prot_motif = "{} {} {} {} {} {}".format(m[0], m[6], m[7], m[8], m[9], m[10])
    motif_inds[prot_motif] = ind

refe_motif = motifs[lineno - 1].split()
print(refe_motif)

refe_pdb = np.load(f"../database/pdbnpy/{refe_motif[0]}.npy")
refe_pdb = refe_pdb[refe_pdb["model"] == 1]
refe_nuc = refe_pdb[refe_pdb["chain"] == refe_motif[1].encode()]
refe_prot = refe_pdb[refe_pdb["chain"] == refe_motif[6].encode()]

desc = calc_stacking_descriptors(refe_motif[7], refe_prot)
desc.update(calc_stacking_descriptors(refe_motif[9], refe_prot))
refe_desc1 = desc[int(refe_motif[8])]
refe_desc2 = desc[int(refe_motif[10])]
refe_desc = np.concatenate((refe_desc1, refe_desc2))
refe_desc_no_ca = np.concatenate((refe_desc1[1:], refe_desc2[1:]))

cbdis = norm(refe_desc[5] - refe_desc[1])
cgdis = norm(refe_desc[6] - refe_desc[2])

found_seqsigs = {}  # seqsig-to-closest-RMSD
unbound_positive_indices = {}
false_positive_seqsigs = {}


def compare_seqsig(seqsig1, seqsig2):
    identities = 0
    for c1, c2 in zip(seqsig1, seqsig2):
        if c1 == c2:
            identities += 1
        elif c1 == "-" or c2 == "-":
            identities += 0.5
    if identities - 3 >= 24:  # 80 %
        return True
    return False


pdbs = glob.glob("../database/pdbnpy/*.npy")  # 2443
for pdbf in sorted(pdbs):
    code = pdbf[pdbf.rindex("/") + 1 : pdbf.rindex(".")]
    # if code not in ("1B7F", "4QQB"):
    #    continue
    # if code != "1B7F":
    #    continue
    # if code not in ("1B7F", "3TRZ", "3TS0", "3TS2", "4QQB", "5UDZ", "5YTT", "6A6J"):
    #    continue
    print(code)
    pdb = np.load(pdbf)
    pdb = pdb[pdb["model"] == 1]
    chains = np.unique(pdb["chain"])
    for chain in chains:
        pdbc = pdb[pdb["chain"] == chain]
        chain_ca = pdbc[pdbc["name"] == b"CA"]
        desc = {}
        for resname in ("PHE", "TYR", "TRP", "HIS"):
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
            rms = rmsd[n]
            if rms > 1:
                break

            resid1, resid2 = resids[p1[n]], resids[p2[n]]
            new_motif = "{} {} {} {} {} {}".format(
                code,
                chain.decode(),
                resnames[resid1].decode(),
                resid1,
                resnames[resid2].decode(),
                resid2,
            )
            if motif_inds.get(new_motif, None) == lineno - 1:
                continue

            new_seqsig = get_seqsig(resid1, resid2, chain_ca)

            hit_seqsigs = []
            if new_seqsig in seqsigs or new_seqsig in found_seqsigs:
                hit_seqsigs.append(new_seqsig)
            else:
                for seqsig in seqsigs:
                    is_close = compare_seqsig(seqsig, new_seqsig)
                    if is_close:
                        hit_seqsigs.append(seqsig)

            if len(hit_seqsigs):
                # seqsig (and/or a close homolog) does bind dinucleotide
                for hit_seqsig in hit_seqsigs:
                    curr_best_rmsd = found_seqsigs.get(hit_seqsig, 9999)
                    if rms < curr_best_rmsd:
                        found_seqsigs[hit_seqsig] = rms

            else:
                # seqsig does not bind dinucleotide

                curr_best_rmsd = false_positive_seqsigs.get(new_seqsig, 9999)
                if rms < curr_best_rmsd:
                    false_positive_seqsigs[new_seqsig] = rms


fp1 = open(truepos_outfile, "w")

fp2 = open(falsepos_outfile, "w")

found_seqsig_list = sorted(found_seqsigs.items(), key=lambda kv: kv[1])

for k, v in found_seqsig_list:
    print(k, "%.5f" % v, file=fp1)

false_positive_seqsig_list = sorted(
    false_positive_seqsigs.items(), key=lambda kv: kv[1]
)

for k, v in false_positive_seqsig_list:
    print(k, "%.5f" % v, file=fp2)
