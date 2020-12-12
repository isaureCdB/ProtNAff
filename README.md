# ProtNAff
Create filters to select Protein - Nucleic acids structures from PDB files, and build libraries of protein-boud RNA fragments 

--------------------------------------------------------------------------
ProtNAff is a pipeline to:
--------------------------------------------------------------------------
1_ Clean-up and parse all NA-protein structures from the PDB into ensembles of small information units in a single file

2_ Search for sets of NA-protein structures with highly cutomisable combinations of criteria

3_ Create RNA/DNA 3D fragment libraries extracted from those sets of structures.

4_ Perform statistics on customised features of such libraries

1 is the necessary first step, 2 and 3 can be done independantly, 4 must be done after 3.
The output of 1 for the PDB at a certain time will soon be downloadable from the LORIA website, allowing to do directly 2 and/or 3(+/-4) directly.

--------------------------------------------------------------------------
Work in progress
--------------------------------------------------------------------------
The main idea of ProtNAff is to provide a highly versatile pipeline to cover as many usages as possible.
This is intended to be a dynamic collaborative work:

_ If you need a specific feature that you can't find or don't know how to add in the current pipeline, please contact us, and we will do our best to include it!

_ If you added some feature that you think can be usefull to others, please fill free to add a new branch in this repository.

--------------------------------------------------------------------------
Steps description
--------------------------------------------------------------------------
--------------------------------------------------------------------------
1. Creation of the NA-protein database
--------------------------------------------------------------------------
The creation of the database is done by:
./create_database.sh [pdbcodes.list] ["dna"/"rna"]

The steps performed by create_database.sh are the following:

_ download all protein-NA structures from the PDB

_ extract relevant information units on each structure (resolution, NA type, etc) and store it into one single easiliy-searchable json file.

_ clean up each structure (add missing atoms, list incomplete nucleotides, list HETATM, etc)

_ characterize the interface (sugar/phosphate/base - protein distances, water contacts, etc)

_ use the 3DNA program [1] for NA structure description that gives exhaustive data in easily parsable Json format

_ rearrange the data per nucleotide (eg. “nucl 5 to 15 make a stem-loop” → “nucl 5 is at position 1 in an 11-nucl stem-loop”)

[1] X-J Lu & WK Olson. 3DNA: a software package for the analysis, rebuilding and visualization of three‐dimensional nucleic acid structures. Nucleic Acids Research (2003) 31(17), 5108-21

The output is a description, at the PDB structure levels, of the full set of protein-bound NAs from the PDB, in a single Json file.


--------------------------------------------------------------------------
2. Requests to select (parts of) structures from the database
--------------------------------------------------------------------------
To use querries you need to use ./create_benchmark/get_benchmark.sh script

The idea is to create an architecture of folders in which you will have for each pdb if in the Json file the result of the query.

The query can ask questions on all informations of the Json file :

_ if nucleotides are in contact with protein
_ if nucleotides have a specific secondary structure
_ etc

You can use exemple of queries which are in ./create_benchmark/queries or construct your own queries.


--------------------------------------------------------------------------
3. Create a fragment library
--------------------------------------------------------------------------
./create_frag_library.sh ["dna"/"rna"]

The protein-bound NA structures are cut into fragments (here trinucleotides), that are pooled by sequence and clustered by pairwise RMSD.

To increase the number of fragments per sequence motif, the full structures are mutated T/U > C and G > A (ex. AGU becomes AAC), and the representative fragments after clustering are mutated back into the 8 possible "de-mutated" sequence (AAC => AAC, GAC, GGC, GGU, AGC, AGU, AAU, GAU).

To accelerate clustering, the fragments are first converted into ATTRACT's coarse-grained representation. The resulting clustering scheme is then applied to the all-atom fragments.


--------------------------------------------------------------------------
4. Requests to compute statistics on a fragments library
--------------------------------------------------------------------------
The fragment library is converted into a flat table for quick search.

Queries on the fragment library use 3 dictionaries:
_ the data dictionnary, provided by structures.json (written when creating the database, at step 1)
_ the chainschema (see make_chainschena.py) that describes the format of the data and how it can be queried
_ the query "variables" dictionnary (see query.py), that contains description of the subdata we are interested in

For details on the data formats, see make_chainschena.py
For an example of query, see the jupyter-notebook exemple_stats.ipynb
