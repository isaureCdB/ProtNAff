import numpy as np
import os
import argparse
p = argparse.ArgumentParser()

p.add_argument("library_pattern", help="""e.g. X=library/CAC-lib. Assumes the existence of:
X.list 
X-conformer.list X-conformer.npy 
X-intracluster.list X-intracluster.npy 
X-secondary.list X-tertiary.list X-quaternary.list X-quasi-unique.list"""
)
p.add_argument("pdb_code", help="PDB code that is to be replaced in the library")
p.add_argument("output", help="Output for the modified library. Contains coordinates unless --write-conformers is set.")
p.add_argument("--max-rmsd", default=1.5, type=float, help="Maximum RMSD to add from the quaternary list. Fragments with a higher RMSD are considered unmodelable.")
p.add_argument("--min-rmsd", default=0.7, type=float, help="Minimum RMSD to consider for replacement from the tertiary list. Fragments with a lower RMSD are considered already good enough without replacement")
p.add_argument("--min-improvement", default=0.2, type=float, help="Minimum RMSD improvement for replacement from the tertiary list.")
p.add_argument("--quasi-unique", default=False, action="store_true", help="Consider quasi-unique replacements for unmodelable fragments. This will make the docking not fully blind")
p.add_argument("--write-conformers", action="store_true", help="Instead of coordinates, write conformer indices (in X-conformer.list, starting from 0) for each modified library structure")

args = p.parse_args()
X = args.library_pattern
pdb_code = args.pdb_code.upper()[:4]

conformer_coors = np.load(X+"-conformer.npy")
intracluster_coors = np.load(X+"-intracluster.npy")

conformers = []
with open(X+"-conformer.list") as f:
    for l in f:
        conformers.append(int(l))

intracluster = []
with open(X+"-intracluster.list") as f:
    for l in f:
        intracluster.append(int(l))

secondary = {}
with open(X+"-secondary.list") as f:
    for l in f:
        ll = l.split()
        if ll[0] != pdb_code:
            continue
        secondary[int(ll[1])] = int(ll[2])

nelim1 = 0
nrepl2 = 0
nrepl3 = 0
nrepl4 = 0
nrepl_quasi = 0

coors = []
conformer_list = []

# Replace primary list with secondary
with open(X+".list") as f:
    for l in f:
        ll = l.split()
        frag = int(ll[0])
        conf = conformers.index(frag)
        if len(ll) == 1 or ll[1] != pdb_code:
            coors.append(conformer_coors[conf])
            conformer_list.append(conf)
        elif secondary.get(frag) is not None:
            frag2 = secondary[frag]
            c = intracluster.index(frag2)
            coors.append(intracluster_coors[c])
            conformer_list.append(conf)
            nrepl2 += 1
            nelim1 += 1
        else:
            nelim1 += 1
if nelim1 > 0:
    print(f"{nelim1} structures eliminated from the primary list")
if nrepl2 > 0:
    print(f"{nrepl2} structures added from the secondary list")

with open(X+"-tertiary.list") as f:
    for l in f:
        ll = l.split()
        if ll[0] != pdb_code:
            continue
        primary_index = int(ll[1])
        new_frag = int(ll[2])
        new_rmsd = float(ll[3])
        old_rmsd = float(ll[4])
        if old_rmsd < args.min_rmsd:
            continue
        if old_rmsd - new_rmsd < args.min_improvement:
            continue
        try:
            conf = intracluster.index(new_frag)
        except ValueError:
            conf = conformers.index(new_frag)
            if conf in conformer_list:
                continue
            coors.append(conformer_coors[conf])
            conformer_list.append(conf)
        else:
            coors.append(intracluster_coors[conf])
        nrepl3 += 1

if nrepl3 > 0:
    print(f"{nrepl3} structures added from the tertiary list")

quasi_unique = []
if args.quasi_unique:
    with open(X+"-quasi-unique.list") as f:
        for l in f:
            ll = l.split()
            if ll[0] != pdb_code:
                continue
            old_frag = int(ll[1])
            new_frag = int(ll[2])
            quasi_unique.append(old_frag)
            try:
                conf = intracluster.index(new_frag)
            except ValueError:
                conf = conformers.index(new_frag)
                if conf in conformer_list:
                    continue
                coors.append(conformer_coors[conf])
                conformer_list.append(conf)
            else:
                coors.append(intracluster_coors[conf])
            nrepl_quasi += 1

with open(X+"-quaternary.list") as f:
    for l in f:
        ll = l.split()
        if ll[0] != pdb_code:
            continue
        old_frag = int(ll[1])
        new_frag = int(ll[2])
        new_rmsd = float(ll[3])
        if old_frag in quasi_unique:
            continue
        if new_rmsd > args.max_rmsd:
            print("Give up on structure: RMSD {:.3f}".format(new_rmsd))
            continue
        try:
            conf = intracluster.index(new_frag)
        except ValueError:
            conf = conformers.index(new_frag)
            if conf in conformer_list:
                continue
            coors.append(conformer_coors[conf])
            conformer_list.append(conf)
        else:
            coors.append(intracluster_coors[conf])
        nrepl4 += 1

if nrepl4 > 0:
    print(f"{nrepl4} structures added from the quaternary list")

if nrepl_quasi > 0:
    print(f"{nrepl_quasi} structures added from the quasi-unique list")

dest = args.output
if args.write_conformers:
    print(f"Write conformer indices (in {X}-conformer.list) to {dest}")
    with open(dest, "w") as f:
        for conf in conformer_list:
            print(conf, file=f)
else:
    if not dest.endswith(".npy"):
        dest += ".npy"
    if not (nelim1 or nrepl2 or nrepl3 or nrepl4 or nrepl_quasi):
        print("No modifications were made")
        src = X + ".npy"
        src2 = os.path.abspath(src)
        if os.path.exists(dest):
            os.remove(dest)
        try:
            os.link(src, dest)
            txt = "Hardlink"
        except OSError:
            os.symlink(src2, dest)
            txt = "Symbolic link"
        print(f"{txt} created: {dest} => {src}")
    else:
        coors = np.stack(coors)
        np.save(dest, coors)
        print(f"Write modified library coordinates to {dest}")

