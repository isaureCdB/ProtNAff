set -u -e

filelist=$1
na=$2

cat /dev/null > parse_pdb_initial.err
for inp in `cat $filelist`; do
  name=${inp%%.pdb}
  outf=$name-iniparse.pdb
  if [ -s $outf ] || [ ! -s $inp ]; then
    echo "skip $name" >> /dev/stderr
    continue
  fi
  chain=`awk '$1=="ATOM"||$1=="HETATM"{print substr($0,22,1);exit}' $inp`
  echo "------------------------------------ "$inp $chain >> /dev/stderr
  echo "------------------------------------ "$inp >> parse_pdb_initial.err
  set +e
  first=`awk '$1=="ATOM"||$1=="HETATM"{printf "%i", substr($0,23,4); exit}' $inp`
  python2 $ATTRACTDIR/../allatom/aareduce.py $inp $outf --$na --heavy \
  --manual --modbase --patch $first None --startres $first --chain $chain \
  --mutate $ATTRACTDIR/../allatom/$na-mutate.list > $name.mapping
  #2>> parse_pdb_initial.err
  if [ -s $outf ]; then
    echo $outf
  fi
  set -e
done
