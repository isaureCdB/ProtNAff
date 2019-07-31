#!/bin/bash

#Author Isaure Chauvot de Beauchene, CNRS, LORIA. 2018

wd=`pwd`
d=`dirname "$0"`
SCRIPTS="$wd/$d/search_frag_library/"

set -u -e

python3 $SCRIPTS/json2npy.py fragments_clust.json fragments_clust.npy

if [ ! -s  /tmp/pdb_chain_pfam.tsv ]; then
    wget ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_pfam.tsv.gz
    gunzip pdb_chain_pfam.tsv.gz
    mv pdb_chain_pfam.tsv /tmp/
python3 $SCRIPTS/add_pfam_uniprot.py

for i in 30 40 50 70 90 95 100; do
    if [ ! -s /tmp/bc-$i.out ];then
        wget ftp://resources.rcsb.org/sequence/clusters/bc-$i.out
        mv bc-$i.out /tmp/
    fi
done
python3 $SCRIPTS/add_seqclust.py
python3 $SCRIPTS/make_chainschema.py
python3 -u $SCRIPTS/build-redundancy-masks.py | tee build-redundancy-masks.stats
