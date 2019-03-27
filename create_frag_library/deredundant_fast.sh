#!/bin/bash
m=$1
dr=$2 # RMSD cutoff for redundancy
a=$3

>&2 echo "Discarde redundant fragments of $m with RMSD < $dr A toward another fragment"

d=`dirname "$0"`

set -u -e

# Remove trinucleotides with inter-nucleotide clashes
# outputs:
#   $m-aa.noclash : list of indices of non-clashing fragments
#   /tmp/$m-aa-noclash.npy : coordinates of non-clashing fragments
if [ ! -s $m-aa.noclash ];then
  echo '#discard clashing fragments $m'
  n=`awk 'BEGIN{j=1}$3=="P"&&NR>1{print NR-j; j=NR}END{print NR-j+1}' templates/$m.pdb`
  $d/discard_clashing_fragments.py $m-all-aa.npy $n 2 /tmp/$m-aa-noclash.npy $m-aa.noclash
fi

# Fit all conformers on the first one
# output:
#   $m-aa-fit.npy : coordinates of fitted fragments
if [ ! -s $m-aa-fit.npy ];then
  echo "#center and fit conformers $m"
  $d/fit_multi_npy.py /tmp/$m-aa-noclash.npy --pdb templates/$m.pdb $m-aa-fit.npy
fi

# Remove redundant conformers (RMSD < ${dr} A)
# output:
#   $m-aa-fit-clust${dr}.npy : coordinates of non-redundant fitted fragments
#   $m-aa-fit-clust${dr} : global indices of non-redundant fragments
#if [ ! -s $m-aa-fit-clust${dr} ];then
  echo "Deredundant the fitted conformers at ${dr}A $m"
  f=$m-aa-fit
  $d/fastcluster_npy.py $f.npy $a 2> fastcluster-$m-$a.log
  $d/concatenate_clusters.py $f-clust$a.0 > $f-clust$a.0-concat
  $d/fastsubcluster_npy.py $f.npy $f-clust$a.0-concat $dr $f-noclash-clust$dr /dev/null 2> fastsubcluster-$m-$dr.log
  $d/select-struct-npy.py $f.npy $f-clust${dr}.npy --structures `awk '{print $4}' $f-noclash-clust${dr}`
  $d/map_cluster.py $f-noclash-clust$dr $m-aa.noclash > $f-clust$dr
  awk '{print $4}' $f-clust${dr} > $m-dr${dr}r.list
#fi
# !!!!! $m-dr${dr}r.list is in absolute indexing, before removing clashes !!!!!!!!!!!!!

# Convert fragments into coarse-grained representation
# output:
#   $m-dr${dr}r.npy : coarse-grained coordinates of non-redundant fitted fragments
if [ ! -s $m-dr${dr}r.npy ];then
  echo "reduce fit-clust${dr}.npy $m"
  $d/reduce-npy.py $f-clust${dr}.npy templates/$m.pdb $m-dr${dr}r.npy > /dev/null
fi

#rm $f-clust${dr}.npy
