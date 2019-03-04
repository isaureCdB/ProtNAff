import numpy as np
import json
fragments = np.load("fragments_clust-aa_missing.npy")
chaindata = json.load(open("chainsmodel_frag_light.json"))
