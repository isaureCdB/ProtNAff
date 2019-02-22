complex=$1
model=$2

name=$complex-$model
inp=/tmp/$name-3dna.pdb
outp1=/tmp/$name-3dna.json
outp2=3dna/$name.snap

x3dna-dssr -i=$inp -o=$outp1 --non-pair --u-turn --po4 --json --idstr=long --more 2> 3dna.errors

#x3dna-snap -i=$inp -o=$outp2 --cutoff=5.0 --get-hbond --t-shape 2>> 3dna.errors
#for j in `ls dssr-*`; do mv $j 3dna/$name-dssr-$j;done

jq . $outp1 > 3dna/$name-dssr.json
echo 3dna/$name-dssr.json >> dssr.list
#echo 3dna/$name.snap >> snap.list

#./scripts/3dna.sh 6AR1 1
exit
and I now would like to run some statistics on fragment conformations, depending on parameters
