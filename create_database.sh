#!/bin/bash

#Author Isaure Chauvot de Beauchene, CNRS, LORIA. 2018

# provide a list of PDB codes of structures to parse
pdbcodes=$1

# provide the type of nucleic acids (dna or rna)
na=$2

wd=`pwd`
d=`dirname "$0"`
SCRIPTS="$d/create_database/"

set -u -e
###############################################################
echo "---------------------------- Download PDBs"
###############################################################
mkdir -p brutPDBs
cd brutPDBs
# Download only files that have not yet been downloaded
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
$SCRIPTS/check_pdb.py brutPDBs/ $pdbcodes corrupted_pdb_files.list \
  tofix.list checked.list splitted.list $SCRIPTS/${na}lib/mutate.list \
   $na
sort -u checked.list > bi; mv -f bi checked.list

##########################################################################
echo "-------------------------- detect NA - protein interface "
##########################################################################
  # Split into chains and models. Creates:
  #   - chainsmodels/xxxxX-y.pdb
  #   - chainsmodels.json
  # options:
  #    --delete : delete from output json the pdbcode entries not in input list
  #    --replace: delete in output json the entries also in input list
  #    defaults: update output json with new entries in input
$SCRIPTS/interface_pdb_contacts.py 5 brutPDBs chainsmodels checked.list \
  $SCRIPTS/${na}lib/mutate.list chainsmodels.json $na \
   > interface_pdb_contacts.log

##########################################################################
echo "---------------------------------  parse initial pdb"
##########################################################################
  ### Remove/rename atoms from modified bases. Creates:
      #   - cleanPDB/xxxxX-y.pdb
      #   - cleanPDB.list
  # Also checks that every chain has " " or "A" in column 17 (alternative conformations)
  # TODO: remove clashing atoms
$SCRIPTS/clean_rna.py 'chainsmodels/' 'cleanPDB/' cleanPDB.list \
  $SCRIPTS/${na}lib/mutate.list chainsmodels.json \
  clean_rna.json $na > clean.err

  ### Applies aareduce to remove non-NA. Creates:
  #    - parse_pdb_initial.errors
  #    - cleanPDB/xxxxX-y-iniparse.pdb
  # Marks missing atoms by XXX coordinates
$SCRIPTS/parse_pdb_initial.sh cleanPDB.list $na >> clean-iniparse.list
sort -u clean-iniparse.list > bi; mv -f bi clean-iniparse.list


# This section is to create a library of mononucleotide, to be used by the
# pdbcompletion.py script to add missing atoms.
# A mono-nucleotide library created in 2018 is already provided, but you might
# want to update/customise it (e.g. change the resolution = the clustering cutoff)
if False;then
    ##########################################################################
    echo "---------------------------- update mononucleotide library"
    ##########################################################################
    $SCRIPTS/build_monolib.py data/${na}lib clean-iniparse.list data/${na}lib

    ut=U; if [ "$na" == "dna" ];then ut=T; fi
    for a in A C G $ut ; do
      $SCRIPTS/cluster_monolib_npy.sh data/${na}lib/$a-fit.npy data/${na}lib/$a.pdb 0.3
      mv data/${na}lib/$a-fit-clust$cut-sel.npy data/${na}lib/$a.npy
    done
    rm data/${na}lib/$a-fit*.npy
fi


##########################################################################
echo "--------------------------------- Fill-up missing atoms "
##########################################################################
  ### Record nucl with missing atoms. Create:
      #    - excise-pdb-missings.err
      #    - excised.list  : List of files checked for missing atoms
      #                      Some chains are discarded
  ### Add entries in clean_rna.json:
      #   - 'missings' : residues with missing atoms (if too many => deleted)
      #   - 'sequence'
$SCRIPTS/excise-pdb-missings.py cleanPDB clean_rna.json excise.json \
  excised.list $na

  ### Fill-up missing atoms
      #    create: _ parse_pdb.errors
      #             _ cleanPDB/xxxxxX-y-iniparse-aa.pdb
      #             _ clean-iniparse-aa.list
      #             _ still_missing.list
      # Add 5'-Phosphate group (5PHO) to 5'-termini.
      # Some 5PHO can still be missing because all conf in nalib clash
      # Those are written as XXX, and should be indicated in a json file
$SCRIPTS/parse_pdb-missings.sh excised.list $na

### retrieve missing (always clashing) 5PHO
$SCRIPTS/excise-still-missing.py cleanPDB still_missing.list \
  clean-iniparse-aa.list excise.json $na
sort -u clean-iniparse-aa.list > bi; mv -f bi clean-iniparse-aa.list

##########################################################################
echo "-------------------------------- apply 3dna"
##########################################################################
mkdir -p 3dna
### concatenate protein and RNA chains for each model###
$SCRIPTS/3dna.py excise.json $SCRIPTS/3dna.sh
# Convert data per structure into data per nucleotide
$SCRIPTS/3dna_parse_json.py excise.json structures.json $na
