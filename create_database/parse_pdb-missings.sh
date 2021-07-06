#TODO: remove --heavy option

filelist=$1 # missings.list
na=$2

d=`dirname "$0"`

set -u -e
cat /dev/null > clean-iniparse-aa.list
cat /dev/null > parse_pdb.errors
cat /dev/null > still_missing.list

for i in `cat  $filelist`; do
  name=${i%%-excise.pdb}
  outf=$name-aa.pdb
  if [ -s $outf ] && [ -s $name-aa.mapping ]; then
      echo $outf >> clean-iniparse-aa.list
      continue
  fi
  if [ ! -s $i ] ; then
    continue
  fi
  resn=`head -1 $i | awk '{print substr($0,23,4)}'`
  chain=`echo $i|awk -F "-" '{print substr($1,length($1),1)}'`

$d/pdbcompletion.py $i $outf --$na --heavy \
--${na}_chain $chain --modbase --mutate $d/../data/${na}lib/mutate.list \
--patch None $resn None --manual residual > $name-aa.mapping 2> pdbcompletion.err

  echo $outf >> clean-iniparse-aa.list

  if [ ! -f $outf ]; then
      touch $outf
  fi
done
