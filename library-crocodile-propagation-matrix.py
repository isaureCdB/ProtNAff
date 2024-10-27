import itertools
import json
import os
import sys
import numpy as np
import argparse

from tqdm import tqdm
from nefertiti.functions.parse_pdb import parse_pdb

def err(msg):
    print(msg, file=sys.stderr)
    exit(1)

parser = argparse.ArgumentParser()
parser.add_argument("pre_sequence", help="trinucleotide sequence pre")
parser.add_argument("post_sequence", help="trinucleotide sequence post")
parser.add_argument("max_overlap_rmsd", type=float)

parser.add_argument(
    "--rotaconformers",
    help="""Rotaconformer index file pattern.
                    Pattern must contain XXX, which will be replaced for each trinucleotide sequence.
                    Each file must contain a list of filenames (or checksums), one per conformer. 
                    Each filename or checksum corresponds to the array of (3, 3) rotation matrices for that conformer
                    """,
    required=True,
)

parser.add_argument(
    "--rotaconformer-directory",
    help="Directory to prepend to the filenames/checksums in the rotaconformer indices",
)

parser.add_argument("--gridspacing", type=float, default=np.sqrt(3)/3)
args = parser.parse_args()

assert len(args.pre_sequence) == 3
assert len(args.post_sequence) == 3
assert args.pre_sequence[-2:] == args.post_sequence[:2]
assert len([nuc for nuc in args.pre_sequence if nuc in ("A", "C")]) == 3
assert len([nuc for nuc in args.post_sequence if nuc in ("A", "C")]) == 3
sequence1, sequence2 = args.pre_sequence, args.post_sequence

gridspacing = args.gridspacing
max_overlap_rmsd = args.max_overlap_rmsd

translation_grids = []
_max_space = int(np.sqrt(max_overlap_rmsd**2/gridspacing**2)) + 1  # +1 because some offsets are initially in the right direction
for n in range(_max_space+1):
    _space = np.arange(-n, n+1)
    translation_grid = gridspacing * np.stack(np.meshgrid(_space, _space, _space, indexing='ij'),axis=-1).reshape(-1, 3)
    translation_grids.append(translation_grid)

library = {}
preatoms = {}
postatoms = {}
for sequence in (sequence1, sequence2):
    template = parse_pdb(open(f"../../templates/{sequence}-template.pdb").read())
    preatoms[sequence] = (template["resid"] >= 2) 
    postatoms[sequence] = (template["resid"] <= 2)
    #print(sequence, preatoms[sequence].sum(), postatoms[sequence].sum())
    library[sequence] = np.load(f"{sequence}-lib-conformer.npy")    
    assert library[sequence].shape[1] == len(template)

rotaconformer_indices = {}
for sequence in (sequence1, sequence2):
    filename = args.rotaconformers.replace("XXX", sequence)
    with open(filename) as fp:
        rotaconformer_index = json.load(fp)
    if not isinstance(rotaconformer_index, list):
        err(
            f"Sequence {sequence}: '{filename}' is not a list of filenames/checksums"
        )
    if len(rotaconformer_index) != len(library[sequence]):
        err(
            f"Sequence {sequence}: There are {len(library[sequence])} conformers but {len(rotaconformer_index)} rotaconformers"
        )
    if args.rotaconformer_directory:
        for fnr, f in list(enumerate(rotaconformer_index)):
            rotaconformer_index[fnr] = os.path.join(
                args.rotaconformer_directory, f
            )
    rotaconformer_indices[sequence] = rotaconformer_index

compatibility_sequence = sequence1 + sequence2[-1:]
compatibility_matrix = np.load(f"{compatibility_sequence}-compatibility-rmsd.npy")

output_dtype = np.dtype(
    [
        ("conformer", np.uint16),
        ("rotamer", np.uint16),
        ("translation", np.float32, 3),
        ("rmsd", np.float32),
    ],
    align=True,
)

assert preatoms[sequence1].sum() == postatoms[sequence2].sum()
lib1 = library[sequence1][:, preatoms[sequence1]]
lib2 = library[sequence2][:, postatoms[sequence2]]

np.random.seed(0)
confs1 = np.array(list(range(len(lib1))))
np.random.shuffle(confs1)
confs1_sample = confs1[:100]
rotamers1_sample = (1000, 1001)

for conf1 in confs1_sample:
    print(conf1+1)
    
    allowed_conf2 = []
    allowed_conf2_translation_grids = []
    for conf2 in range(len(lib2)):
        overlap_rmsd = compatibility_matrix[conf1, conf2]
        if overlap_rmsd >= max_overlap_rmsd:
            continue
        
        # DO NOT FILTER HERE: filter post-hoc
        #if overlap_rmsd >= 0.5:  # 70 % of the cases
        #    continue
        #if overlap_rmsd >= 0.75:  # 95 % of the cases
        #    continue

        allowed_conf2.append(conf2)
        
        residual_error = max_overlap_rmsd**2 - overlap_rmsd**2
        min_allowed_grid_steps = int(np.sqrt(residual_error / (gridspacing**2))) + 1
        #print(overlap_rmsd, residual_error/(gridspacing**2), min_allowed_grid_steps)
        allowed_translation_grid = translation_grids[min_allowed_grid_steps]
        allowed_conf2_translation_grids.append(allowed_translation_grid)
    if not len(allowed_conf2): continue

    coor1 = lib1[conf1]
    roco_matrices1 = np.load(rotaconformer_indices[sequence1][conf1])
    
    for rotamer1 in rotamers1_sample:
        assert rotamer1 < len(roco_matrices1), (rotamer1, roco_matrices1)

    roco_coors1 = np.einsum("kj,ijl->ikl", coor1, roco_matrices1)
    coms1 = roco_coors1.mean(axis=1)
    
    for rotamer1 in rotamers1_sample:
        outputfile = f"propagation/{sequence1}-{conf1+1}-{rotamer1+1}-{sequence2}.npy"
        results = np.empty(1000000, dtype=output_dtype)
        nresults = 0
        roco_coor1, com1 = roco_coors1[rotamer1], coms1[rotamer1]

        # 1.3 A max_overlap_rmsd:
        # 5 minutes without conformer overlap rmsd filtering (CURRENT VERSION)
        # 1-2 minutes with conformer overlap rmsd filtering, threshold 0.75 A
        # max 20 secs (often 3 secs) with conformer overlap rmsd filtering, threshold 0.5 A

        for conf2, trans_grid in tqdm(zip(allowed_conf2, allowed_conf2_translation_grids), total=len(allowed_conf2)):
            coor2 = lib2[conf2]
            roco_matrices2 = np.load(rotaconformer_indices[sequence2][conf2])
            # to optimize: most conf2 rotamers will be useless!
            roco_coors2 = np.einsum("kj,ijl->ikl", coor2, roco_matrices2)
            offsets = com1[None, :] - roco_coors2.mean(axis=1)                
            offsets_disc = np.round(offsets / gridspacing) * gridspacing
            dif = roco_coors2 + offsets_disc[:, None] - roco_coor1
            rmsds = np.sqrt((dif * dif).sum(axis=2).mean(axis=1))
            mask = (rmsds < max_overlap_rmsd)
            mask_ind = np.where(mask)[0]
            #/to optimize
            if not len(mask_ind):
                continue

            offsets_disc_masked = offsets_disc[mask]
            base_dif = roco_coors2[mask] + offsets_disc_masked[:, None] - roco_coor1
            grid_dif = base_dif[:, None, :, :] + trans_grid[None, :, None, :]
            grid_rmsds = np.sqrt((grid_dif * grid_dif).sum(axis=3).mean(axis=2))
            grid_mask = (grid_rmsds < max_overlap_rmsd)
            
            grid_rota_mask, grid_trans_mask  = np.where(grid_mask)
            ncurr_results = len(grid_rota_mask)
            if not ncurr_results:
                continue
            curr_results = results[nresults:nresults+ncurr_results]
            curr_results["conformer"] = conf2
            curr_results["rotamer"] = mask_ind[grid_rota_mask]
            curr_results["translation"] = offsets_disc_masked[grid_rota_mask] + trans_grid[grid_trans_mask]
            curr_results["rmsd"] = grid_rmsds[grid_mask]
            nresults += ncurr_results

        np.save(outputfile, results[:nresults])
