# launch from the main ProtNAff directory. 
# Takes about 20 hours on 8 cores
for motif in AAA AAC ACA ACC CAA CAC CCA CCC; do
    for pos in 0 1; do
        for dihexmod0 in 1 9 17 25; do
            echo $motif $pos $dihexmod0
            last=$((dihexmod0+7))
            for dihexmod in $(seq $dihexmod0 $last); do
                python dinucleotide-completeness/dinucleotide-completeness.py $motif $pos $dihexmod > TMP-$motif-$pos-$dihexmod &
            done
            wait
        done
        cat TMP-$motif-$pos-* | sort -nk3 > database/dinucleotide-completeness-$motif-$pos.txt
        rm -f TMP-$motif-$pos
    done
done