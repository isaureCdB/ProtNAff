#!/bin/bash

name=${1%%.npy}   # /data/${na}lib/A-fit.npy
template=$2       # /data/${na}lib/A.pdb
cut=$3

d=`dirname "$0"`

# fast clustering based on 3 atoms
$d/select-struct-npy.py ${name}.npy ${name}-3at.npy \
  --atnames "P" "O3'" "C5" --template $template

cut1=$((cut+1))
echo $cut $cut1
echo "$d/fastcluster_npy.py ${name}-3at.npy $cut1 --chunk 200"
$d/fastcluster_npy.py ${name}-3at.npy $cut1 --chunk 200 #${name}-clust3.npy

# Slow sub-clusteringi nall-atoms
$d/subcluster-npy.py ${name}.npy $name-3at-clust$cut1 $cut

awk 'NF>12{print $4}' $name-clust$cut > $name-clust$cut-sel
$d/select-struct-npy.py ${name}.npy $name-clust$cut-sel.npy --structures `cat $name-clust$cut-sel`
