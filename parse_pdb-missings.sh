#TODO: remove --heavy option
filelist=$1 # missings.list
na=$2
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
    echo "$i does not exist "
    continue
  fi
  echo $i $outf
  resn=`head -1 $i | awk '{print substr($0,23,4)}'`
  chain=`echo $i|awk -F "-" '{print substr($1,length($1),1)}'`

  python2 $ATTRACTDIR/../allatom/aareduce.py $i $outf --$na --heavy --nalib \
  --chain $chain --modbase --mutate $ATTRACTDIR/../allatom/$na-mutate.list \
  --patch $resn None --manual > $name-aa.mapping

#  a=`grep XXX $outf|wc -l`
#  if [ "$a" -gt 0 ]; then
#    mv $outf $name-aa-miss5PHO.pdb
#    echo $name-aa-miss5PHO.pdb >> still_missing.list
#  else
  echo $outf >> clean-iniparse-aa.list
#  fi

  if [ ! -f $outf ]; then
      touch $outf
      echo "missing $outf"
  fi
done

#5I4AD
