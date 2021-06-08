name=$1 # $m-dr${dr}r
m=$2  # sequence motif
cut1=$3 # Tight clustering cutoff (1.0)
cut2=$4 # Large clustering cutoff (2.0)

d=`dirname "$0"`

source $MY_CONDA/etc/profile.d/conda.sh
conda activate protnaff
set -u -e

# Cluster fragments by RMSD.
a=`echo $cut1 + 1|bc`
>&2 echo "-------------------------------------------------"
>&2 echo "Cluster $m large ($a A); then tight ($cut1 A)"
>&2 echo "-------------------------------------------------"

# Hierarchical clustering (for speed & memory):
#   first cluster at a = (cut + 1) A,
>&2 echo "cluster $m large $m ($a A)"
python3 $d/fastcluster_npy.py $name.npy $a --chunk 200
# concatenate small clusters in one, to diminish edge effect
python3 $d/concatenate_clusters.py $name-clust$a > $name-clust$a-concat

#   then sub-cluster each ${a}A-cluster members at $cut A
>&2 echo "cluster $m tight $m ($cut1 A)"
python3 $d/subcluster-npy.py $name.npy $name-clust${a}-concat $cut1 > /dev/null

# Extract coordinates of cluster-centers
awk '{print $4}' $name-clust$cut1 > $name-clust$cut1.centers
$d/select-struct-npy.py $name.npy $name-clust$cut1.npy --structure `cat $name-clust$cut1.centers`

# !!!!! $name.list is in absolute indexing, before removing clashes !!!!!
$d/select-lines.py $name.list $name-clust$cut1.centers > $name-clust$cut1.list

# get indices of cluster centers in $m-aa.noclash
$d/select-lines.py $m-aa.noclash $name-clust$cut1.list --indices --lines > $name-clust$cut1-noclash.list

$d/select-struct-npy.py $m-aa-fit.npy $name-clust$cut1-aa.npy  --structure `cat $name-clust$cut1-noclash.list`

# Super-cluster the cluster-centers
python3 $d/cluster-npy.py $name-clust$cut1.npy $cut2 > /dev/null

rm $name-clust$a $name-clust$a.npy $name-superclust${cut1} $name-clust$cut1.centers
