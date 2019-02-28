m=$1
dr=$2
cut=$3

d=`dirname "$0"`

echo "Convert $m library into coarse-grained"

sed 's/conf/confr/' $m-clust${cut}-aa.list > $m-clust${cut}\r.list
python $d/reduce.py $m-clust${cut}-aa.list $m-clust${cut}\r.list --batch --rna
$d/pdb2npy.py $m-clust${cut}\r.list --outp $m-clust${cut}r.npy --list
