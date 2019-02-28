# provide a list of PDB codes of structures to parse
pdbcodes=$1

# provide the type of nucleic acids (dna or rna)
na=$2

wd=`pwd`
d=`dirname "$0"`
NAFRAG="$d/scripts"
DATA="$d/data/"

set -u -e
#if false;then
###############################################################
echo "---------------------------- Download PDBs"
###############################################################
if [ ! -d brutPDBs ];then mkdir brutPDBs; fi
cd brutPDBs
for i in `cat ../$pdbcodes|awk '{print toupper($0)}'`; do
    if [ ! -s $i.pdb ] && [ ! -s $i.pdb.bz2 ] ;then
      echo "downloading $i"
      pdb_download $i > /dev/null 2> /dev/null
      sed -i 's/SE   MSE/ SD  MSE/' $i.pdb
      sed -i 's/MSE/MET/' $i.pdb
    elif [ -s $i.pdb.bz2 ];then
      if [ 1 -gt $(ls  ./cleanPDB/${i}[A-Z]-1.pdb 2>/dev/null | wc -w) ]; then
        bunzip2 $i.pdb.bz2
      fi
    fi
done

cd $wd
mkdir -p chainsmodels/
mkdir -p cleanPDB/
mkdir -p interface/
mkdir -p 3dna/
##########################################################################
echo "-------------------------- check pdb"
##########################################################################
$NAFRAG/check_pdb.py brutPDBs/ $pdbcodes corrupted_pdb_files.list \
  tofix.list checked.list splitted.list $DATA/${na}lib/mutate.list \
   $na > log
sort -u checked.list > bi; mv -f bi checked.list

##########################################################################
echo "-------------------------- detect NA - protein interface "
##########################################################################
  # Split into chains and models. Creates:
  #   - chainsmodels/xxxxX-y.pdb
  #   - chainsmodels.json
  #!!! It deletes pdb if already present in chainsmodels.json !!
$NAFRAG/interface_pdb_contacts.py 5 brutPDBs chainsmodels checked.list \
  $DATA/${na}lib/mutate.list chainsmodels.json $na \
   > interface_pdb_contacts.log

##########################################################################
echo "---------------------------------  parse initial pdb"
##########################################################################
  ### Remove/rename atoms from modified bases. Creates:
      #   - cleanPDB/xxxxX-y.pdb
      #   - cleanPDB.list
  # Also checks that every chain has " " or "A" in column 17 (alternative conformations)
  # TODO: remove clashing atoms
$NAFRAG/clean_rna.py 'chainsmodels/' 'cleanPDB/' cleanPDB.list \
  $DATA/${na}lib/mutate.list chainsmodels.json \
  clean_rna.json $na > clean.err

  ### Applies aareduce to remove non-NA. Creates:
  #    - parse_pdb_initial.errors
  #    - cleanPDB/xxxxX-y-iniparse.pdb
  # Marks missing atoms by XXX coordinates
$NAFRAG/parse_pdb_initial.sh cleanPDB.list $na >> clean-iniparse.list
sort -u clean-iniparse.list > bi; mv -f bi clean-iniparse.list

# if False;then

# This section is to create a library of mono-nucleotide, to be used by the
# pdbcompletion.py script to add missing atoms.
# A mono-nucleotide library created in 2018 is already provided, but you might
# want to update/customise it (e.g. change the resolution = the clustering cutoff)
##########################################################################
echo "---------------------------- build mononucleotide library"
##########################################################################
$NAFRAG/build_monolib.py templates clean-iniparse.list monolib

rd=R
ut=U
if [ "$na" == "dna" ];then rd=D; ut=T; fi

for a in A C G $ut ; do
  cluster_monolib.sh monolib/$rd$a-fit.npy 0.3 templates/$rd$a.pdb
done
#fi

##########################################################################
echo "--------------------------------- Fill-up missing atoms "
##########################################################################
  ### Records nucl with missing atoms. Creates
      #    - excise-pdb-missings.err
      #    - excised.list  : List of files checked for missing atoms
      #                      Some chains are discarded
  ### Adds entries in chainsmodels.json:
      #   - 'missings' : residues with missing atoms (if too many => deleted)
      #   - 'sequence'
$NAFRAG/excise-pdb-missings.py cleanPDB clean_rna.json excise.json \
  excised.list $na

  ### Fill-up missing atoms
      #    creates: _ parse_pdb.errors
      #             _ cleanPDB/xxxxxX-y-iniparse-aa.pdb
      #             _ clean-iniparse-aa.list
      #             _ still_missing.list
      # Add 5'-Phosphate group (5PHO) to 5'-termini.
      # Some 5PHO can still be missing because all conf in nalib clash
      # Those are written as XXX, and should be indicated in a json file
$NAFRAG/parse_pdb-missings.sh excised.list $na

### retrieve missing (always clashing) 5PHO
$NAFRAG/excise-still-missing.py cleanPDB still_missing.list \
  clean-iniparse-aa.list excise.json $na
sort -u clean-iniparse-aa.list > bi; mv -f bi clean-iniparse-aa.list

##########################################################################
echo "-------------------------------- apply 3dna"
##########################################################################
### concatenate protein and RNA chains for each model###
mkdir -p 3dna
$NAFRAG/3dna.py excise.json $NAFRAG/3dna.sh
$NAFRAG/3dna_parse_json.py excise.json x3dna.json $na

##########################################################################
echo "-------------------------------- parse pdb files"
##########################################################################
### Mutate G -> A and U/T -> C to increase the number of conformers per purine-pyrimidine motif
$NAFRAG/parse_pdb_AC.sh clean-iniparse-aa.list $na >> pdbfiles-AC.list
sort -u pdbfiles-AC.list > bi; mv -f bi pdbfiles-AC.list

##########################################################################
echo "-------------------------------- cut into fragments"
##########################################################################
mkdir -p PDBs
mkdir -p trilib
$NAFRAG/fragmt-from-AC.py x3dna.json fragments_ori.json structures.json $na motifs.list 'cleanPDB'

cd trilib
mkdir -p templates/
ln -sf ../fragments_ori.json
ln -s ../motifs.list

d=`dirname "$0"`
$NAFRAG/create_templates.py templates $na

##########################################################################
echo "-------------------------------- fragments clustering"
##########################################################################
dr=0.2 # redundancy cutoff
c1=1.0 # tight clustering cutoff
c2=2.0 # large clustering cutoff

# Deredundant fragments at $dr A RMSD
for m in `cat motifs.list`; do
  $NAFRAG/deredundant.sh $m $dr > deredundant-$m.log &
done
wait

# Cluster fragments at $c1 A RMSD
# Clustering the $c1\A-cluster centers at $c2 A
for m in `cat motifs.list`; do
  $NAFRAG/clusterfrag_npy.sh $m-dr${dr}r $m $c1 $c2 > clusterfrag_npy_tight-$m.log &
done
wait

for m in `cat motifs.list`; do
  $NAFRAG/create_libraries.sh $m $dr $c1 # Write PDB for each cluster center
  $NAFRAG/reduce_libraries.sh $m $dr $c1 # Convert to coarse-grained representation
done

# Assign each fragment to its clusters in the json file
$NAFRAG/assign_clusters.py fragments.json $na \
  --clustfiles "aa-fit-clust$dr" "dr0.2r-clust$c1" "dr0.2r-clust$c1-clust$c2" \
  --clustnames "clust$dr" "clust$c1" "clust$c2"

##########################################################################
echo "-------------------------------- mutate back into all sequences"
##########################################################################
$NAFRAG/mutate-AC-libraries_npy.py $na "dr${dr}r-clust$c1"
