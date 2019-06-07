#!/usr/bin/env python3
import numpy as np
import json, sys
fragments = np.load(sys.argv[1]) #fragments_clust-aa_missing.npy
chaindata = json.load(open(sys.argv[2])) #chainsmodel_frag.json
