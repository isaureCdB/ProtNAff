# Creates fragments_demut.json ( frag-index : [frag-index-old, sequence, [mutation status from ori]] )

d=`dirname "$0"`

mutated=$1 #mutated.list
c1=$2

for motif in `cat $mutated`; do
    echo $motif
    motiflist=$motif-clust$c1.list
    aalist=$motif-aa-clust$c1-aa.list
    reducedlist=$motif-clust$c1\r.list
    sed "s/\.pdb/-aa\.pdb/g" $motiflist > $aalist
    sed "s/\.pdb/r\.pdb/g" $motiflist > $reducedlist
    python $d/pdbcompletion.py $motiflist $aalist --batch --rna --heavy --nalib --termini --manual
    a=`cat $motiflist |wc -l`
    for i in `seq $a`;do
        k=$RANDOM
        conf=`awk -v i=$i 'NR==i' $aalist`
        egrep -v XXXX $conf > /tmp/tmp-$k; mv -f /tmp/tmp-$k $conf
        done
    python $ATTRACTTOOLS/reduce.py $aalist $reducedlist --batch --rna
done
