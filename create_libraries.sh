m=$1
dr=$2
cut1=$3

name=$m-dr${dr}r-clust$cut1

d=`dirname "$0"`

echo "create PDB fragments (for mutating) from clust1.0, $m"

$d/select-struct-npy.py $m-aa-fit.npy  $name-aa.npy  --struct `awk '{print $4}' $name`

$d/npy2pdb.py $name-aa.npy  templates/$m.pdb   >  /tmp/$name-aa.pdb

mkdir -p $m
$d/splitmodel.py /tmp/$name-aa.pdb  $m/conf-aa > $m-clust${cut1}-aa.list
