motif=$1
for i in $(seq 32); do
    python3 trinucleotide-completeness/triinucleotide-completeness.py $motif $i >> database/trinucleotide-completeness-$motif-$i
    echo $motif $i
done