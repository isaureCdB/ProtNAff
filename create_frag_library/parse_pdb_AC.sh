#TODO: remove --heavy option when hydrogens added to monomucl library

filelist=$1 #clean-iniparse-aa.list
d=`dirname $0`
na=$2
mutationlist=$3

cat /dev/null > parse_pdb_AC.log
for i in `cat $filelist`; do
  outf=${i%%.pdb}-AC.pdb
  if [ -s $outf ]; then continue; fi
  if [ -s $i ] ; then
    echo $outf
    echo $outf >> parse_pdb_AC.log
    resn=`head -1 $i | awk '{print substr($0,23,4)}'`

    $d/pdbcompletion.py $i $outf --$na --heavy --mutate $mutationlist \
    --patch None $resn None > /dev/null 2>> parse_pdb_AC.log

  else
    echo  "$i does not exist"
    touch $outf
  fi
done
