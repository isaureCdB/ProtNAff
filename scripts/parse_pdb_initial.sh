set -u -e

filelist=$1
na=$2

d=`dirname "$0"`

cat /dev/null > parse_pdb_initial.err
for inp in `cat $filelist`; do
  name=${inp%%.pdb}
  outf=$name-iniparse.pdb
  if [ -s $outf ] || [ ! -s $inp ]; then
    #echo "skip $name" >> /dev/stderr
    continue
  fi
  chain=`awk '$1=="ATOM"||$1=="HETATM"{print substr($0,22,1);exit}' $inp`
  echo "------------------------------------ "$inp $chain >> /dev/stderr
  echo "------------------------------------ "$inp >> parse_pdb_initial.err
  first=`awk '$1=="ATOM"||$1=="HETATM"{printf "%i", substr($0,23,4); exit}' $inp`

  $d/pdbcompletion.py $inp $outf --$na --heavy \
  --manual --modbase --patch None $first None --startres $first \
  --mutate $d/../data/${na}lib/mutate.list > $name.mapping

  if [ -s $outf ]; then
    echo $outf
  fi
done
