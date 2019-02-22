name=$1 # $m-dr${dr}r
m=$2  # sequence motif
cut1=$3 # Tight clustering cutoff
cut2=$4 # Large clustering cutoff

d=`dirname "$0"`

set -u -e

# Cluster fragments by RMSD.
a=`echo $cut1 + 2|bc`
>&2 echo "-------------------------------------------------"
>&2 echo "Cluster $m large ($a A); then tight ($cut1 A)"
>&2 echo "-------------------------------------------------"

# Hierarchical clustering (for speed & memory):
#   first cluster at a = (cut + 2) A,
>&2 echo "cluster $m large $m ($a A)"
python3 $d/fastcluster_npy.py $name.npy $a --chunk 200

#   then sub-cluster each ${a}A-cluster members at $cut A
>&2 echo "cluster $m tight $m ($cut1 A)"
python3 $d/subcluster-npy.py $name.npy $name-clust${a} $cut1 > /dev/null

# Extract coordinates of cluster-centers
awk '{print $4}' $name-clust$cut1 > $name-clust$cut1.centers
$d/select-struct-npy.py $name.npy $name-clust$cut1.npy --structure `cat $name-clust$cut1.centers`

# !!!!! $name.list is in absolute indexing, before removing clashes !!!!!!!!!!!!!
$d/select-lines.py $name.list $name-clust$cut1.centers > $name-clust$cut1.list

$d/select-struct-npy.py $m-aa-fit.npy $name-clust$cut1-aa.npy  --structure `cat $name-clust$cut1.list`

# get gobal indices of cluster centers

# Super-cluster the cluster-centers
python3 $d/cluster-npy.py $name-clust$cut1.npy $cut2 > /dev/null

#rm $name-clust$a $name-clust$a.npy $name-superclust${cut1} $name-clust$cut1.centers
