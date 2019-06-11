#!/usr/bin/env python3
import numpy as np
import json, sys
fragments = np.load(sys.argv[1]) #fragments_clust.npy
chaindata = json.load(open(sys.argv[2])) #structures.json
