#TODO: remove --heavy option
#set -u -e

filelist=$1 #clean-iniparse-aa.list
d=`dirname $0`
na=$2
rm -f parse_pdb_AC.errors

cat /dev/null > parse_pdb_AC.errors
for i in `cat $filelist`; do
  outf=${i%%.pdb}-AC.pdb
  if [ -f $outf ]; then continue; fi
  if [ -s $i ] ; then
    echo $outf >> parse_pdb_AC.errors
    resn=`head -1 $i | awk '{print substr($0,23,4)}'`

    $d/pdbcompletion.py $i $outf --$na --heavy \
    --mutate $d/../data/${na}lib/mutate.list \
    --patch None $resn None > /dev/null 2> /dev/null

  else
    touch $outf
  fi
done
