m=$1
dr=$2 # RMSD cutoff for redundancy

>&2 echo "Discarde redundant fragments of $m with RMSD < $dr A toward another fragment"

d=`dirname "$0"`

set -u -e

# Remove trinucleotides with inter-nucleotide clashes
# outputs:
#   $m-aa.noclash : list of indices of non-clashing fragments
#   /tmp/$m-aa-noclash.npy : coordinates of non-clashing fragments
if [ ! -s $m-aa-noclash.npy ];then
  echo '#discard clashing fragments $m'
  n=`awk 'BEGIN{j=1}$3=="P"&&NR>1{print NR-j; j=NR}END{print NR-j+1}' templates/$m.pdb`
  $d/discard_clashing_fragments.py ../$m-all-aa.npy $n 2 /tmp/$m-aa-noclash.npy $m-aa.noclash
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
if [ ! -s $m-aa-fit-clust${dr} ];then
  echo "Deredundant the fitted conformers at ${dr}A $m"
  $d/fastcluster_npy.py $m-aa-fit.npy ${dr} --chunk 200 --mapping $m-aa.noclash 2> fastcluster-$m-${dr}.log
  awk '{print $4}' $m-aa-fit-clust${dr} > $m-dr${dr}r.list
fi
# !!!!! $m-dr${dr}r.list is in absolute indexing, before removing clashes !!!!!!!!!!!!!


# Convert fragments into coarse-grained representation
# output:
#   $m-dr${dr}r.npy : coarse-grained coordinates of non-redundant fitted fragments
if [ ! -s $m-dr${dr}r.npy ];then
  echo "reduce fit-clust${dr}.npy $m"
  python $ATTRACTTOOLS/reduce-npy.py $m-aa-fit-clust${dr}.npy templates/$m.pdb $m-dr${dr}r.npy > /dev/null
fi

rm $m-aa-fit-clust${dr}.npy
