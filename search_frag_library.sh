#!/bin/bash

#Author Isaure Chauvot de Beauchene, CNRS, LORIA. 2018

wd=`pwd`
d=`dirname "$0"`
SCRIPTS="$wd/$d/search_frag_library/"

set -u -e

python3 $SCRIPTS/json2npy.py fragments_clust.json fragments_clust.npy
python3 $SCRIPTS/add_pfam_uniprot.py
python3 $SCRIPTS/add_seqclust.py
python3 $SCRIPTS/make_chainschema.py
python3 -u $SCRIPTS/build-redundancy-masks.py | tee build-redundancy-masks.stats