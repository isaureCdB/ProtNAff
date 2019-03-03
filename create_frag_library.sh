#!/bin/bash

#Author Isaure Chauvot de Beauchene, CNRS, LORIA. 2018

# provide the type of nucleic acids (dna or rna)
na=$1

wd=`pwd`
d=`dirname "$0"`
SCRIPTS="$d/create_frag_library/"

##########################################################################
echo "-------------------------------- mutate pdb files"
##########################################################################
### Mutate G -> A and U/T -> C to increase the number of conformers per purine-pyrimidine motif
$SCRIPTS/parse_pdb_AC.sh clean-iniparse-aa.list $na >> pdbfiles-AC.list
sort -u pdbfiles-AC.list > bi; mv -f bi pdbfiles-AC.list

##########################################################################
echo "-------------------------------- cut into fragments"
##########################################################################
mkdir -p PDBs
mkdir -p trilib
$SCRIPTS/fragmt-from-AC.py structures.json fragments.json structures_AC.json $na motifs.list 'cleanPDB'

cd trilib
mkdir -p templates/
ln -sf ../fragments.json
ln -s ../motifs.list

d=`dirname "$0"`
$SCRIPTS/create_templates.py templates $na

##########################################################################
echo "-------------------------------- fragments clustering"
##########################################################################
dr=0.2 # redundancy cutoff
c1=1.0 # tight clustering cutoff
c2=2.0 # large clustering cutoff

# Deredundant fragments at $dr A RMSD
for m in `cat motifs.list`; do
  $SCRIPTS/deredundant.sh $m $dr > deredundant-$m.log &
done
wait

# Cluster fragments at $c1 A RMSD
# Clustering the $c1\A-cluster centers at $c2 A
for m in `cat motifs.list`; do
  $SCRIPTS/clusterfrag_npy.sh $m-dr${dr}r $m $c1 $c2 > clusterfrag_npy_tight-$m.log &
done
wait

for m in `cat motifs.list`; do
  $SCRIPTS/create_libraries.sh $m $dr $c1 # Write PDB for each cluster center
  $SCRIPTS/reduce_libraries.sh $m $dr $c1 # Convert to coarse-grained representation
done

# Assign each fragment to its clusters in the json file
$SCRIPTS/assign_clusters.py fragments.json $na \
  --clustfiles "aa-fit-clust$dr" "dr0.2r-clust$c1" "dr0.2r-clust$c1-clust$c2" \
  --clustnames "clust$dr" "clust$c1" "clust$c2"

##########################################################################
echo "-------------------------------- mutate back into all sequences"
##########################################################################
$SCRIPTS/mutate-AC-libraries_npy.py $na "dr${dr}r-clust$c1"
