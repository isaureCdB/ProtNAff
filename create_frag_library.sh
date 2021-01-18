#!/bin/bash

#Author Isaure Chauvot de Beauchene, CNRS, LORIA. 2018

# provide the type of nucleic acids (dna or rna)
na=$1

wd=`pwd`
x=$PROTNAFF
SCRIPTS="$x/create_frag_library/"

rm -rf templates
set -u -e

if [ ! -s fragments.json ];then
  ##########################################################################
  echo "-------------------------------- mutate pdb files"
  ##########################################################################
  ### Mutate G -> A and U/T -> C to increase the number of conformers per purine-pyrimidine motif
  $SCRIPTS/parse_pdb_AC.sh clean-iniparse-aa.list $na $x/data/mutate-AC.list

  ##########################################################################
  echo "-------------------------------- cut into fragments"
  ##########################################################################
  python $SCRIPTS/fragmt-from-AC.py structures.json fragments.json $na motifs.list 'cleanPDB'
fi

mkdir -p trilib
cd trilib
\cp ../motifs.list .
for m in `cat motifs.list`; do
    if [ ! -f $m-all-aa.npy ];then
        ln -s ../$m-all-aa.npy
    fi
done

echo "create_template"
mkdir -p templates
python $SCRIPTS/create_templates.py $SCRIPTS/../data/${na}lib templates $na

##########################################################################
echo "-------------------------------- fragments clustering"
##########################################################################
dr=0.2 # redundancy cutoff
c1=1.0 # tight clustering cutoff
c2=3.0 # large clustering cutoff

# Deredundant fragments at $dr A RMSD
if [ ! -s  AAA-dr0.2r.npy ];then
  for m in `cat motifs.list`; do
    $SCRIPTS/deredundant_fast.sh $m $dr $c2 > deredundant-$m.log &
  done
  wait
fi

# Cluster fragments at $c1 A RMSD
# Clustering the ${c1}A-cluster centers at $c2 A
if [ ! -s AAA-dr0.2r-clust1.0.npy ] ;then
  for m in `cat motifs.list`; do
    $SCRIPTS/clusterfrag_npy.sh $m-dr${dr}r $m $c1 $c2 > clusterfrag_npy_tight-$m.log &
  done
  wait
fi

for m in `cat motifs.list`; do
  $SCRIPTS/create_libraries.sh $m $dr $c1 # Write PDB for each cluster center
  $SCRIPTS/reduce_libraries.sh $m $dr $c1 # Convert to coarse-grained representation
done

# Assign each fragment to its clusters in the json file
python $SCRIPTS/assign_clusters.py ../fragments.json ../fragments_clust.json $na \
  --clustfiles "aa-fit-clust$dr" "dr0.2r-clust$c1" "dr0.2r-clust$c1-clust$c2" \
  --clustnames "clust$dr" "clust$c1" "clust$c2"

##########################################################################
echo "-------------------------------- mutate back into all sequences"
##########################################################################
# Example :
# You have to be located on the dna/rna lib folder
# with the 8 bases 3-mers (GGG, GGT, GTG, GTT, TGG, TGT, TTG, TTT)
# dr0.2r-clust1.0 is a reference to the xxx-dr0.2r-clust1.0.npy files
python $SCRIPTS/mutate-AC-libraries_npy.py $na "dr${dr}r-clust$c1"
