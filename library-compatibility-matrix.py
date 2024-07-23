from tqdm import tqdm
import itertools

import numpy as np
from nefertiti.functions.parse_pdb import parse_pdb
from nefertiti.functions.superimpose import superimpose_array

sequences = ["".join(seq) for seq in itertools.product(('A','C'), repeat=3)]
library = {}
preatoms = {}
postatoms = {}
for sequence in sequences:
    template = parse_pdb(open(f"../../templates/{sequence}-template.pdb").read())
    preatoms[sequence] = (template["resid"] >= 2) 
    postatoms[sequence] = (template["resid"] <= 2)
    #print(sequence, preatoms[sequence].sum(), postatoms[sequence].sum())
    library[sequence] = np.load(f"{sequence}-lib-conformer.npy")
    assert library[sequence].shape[1] == len(template)

compatibility_sequences = ["".join(seq) for seq in itertools.product(('A','C'), repeat=4)]
for compatibility_sequence in compatibility_sequences:
    sequence1 = compatibility_sequence[:3]
    sequence2 = compatibility_sequence[-3:]
    assert preatoms[sequence1].sum() == postatoms[sequence2].sum()
    lib1 = library[sequence1][:, preatoms[sequence1]]
    lib2 = library[sequence2][:, postatoms[sequence2]]
    overlap_rmsd = np.empty((len(lib1), len(lib2)))
    for conf2nr, conf2 in tqdm(enumerate(lib2), desc=compatibility_sequence, total=len(lib2)):
        _, curr_overlap_rmsds = superimpose_array(lib1, conf2)
        overlap_rmsd[:, conf2nr] = curr_overlap_rmsds
    np.save(f"{compatibility_sequence}-compatibility-rmsd.npy", overlap_rmsd)
