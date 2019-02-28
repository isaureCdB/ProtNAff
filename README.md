# NAfragDB
Multi-Purpose Structural Database of Nucleic-Acid/Protein Complexes

--------------------------------------------------------------------------
NAfragDB is a pipeline to:
--------------------------------------------------------------------------
1_ select 3D structures of NA-protein complexes in the PDB with highly cutomisable combinations of criteria

2_ create RNA/DNA 3D structural libraries extracted from those structures from PDB.

3_ perform statistics on customised features of that library

--------------------------------------------------------------------------
Work in progress
--------------------------------------------------------------------------
The main idea of NAfragDB is to provide a highly versatile pipeline to cover as many usages as possible.
This is intended to be a dynamic collaborative work:

_ If you need a specific feature that you can't find or don't know how to add in the current pipeline, please contact us, and we will do our best to include it!

_ If you added some feature that you think can be usefull to other, please fill free to add a new version in this repository.

--------------------------------------------------------------------------
Creation of the library
--------------------------------------------------------------------------
The creation of the 3D library is done by:
./nalib.sh [pdbcodes.list] ["dna"/"rna"]

The steps of ./nalib are to:

_ download all protein-NA structures from the PDB

_ extract relevant information units on each structure (resolution, NA type, etc) and store it into one single easiliy-searchable json file.

_ clean up each structure (add missing atoms, list incomplete nucleotides, list HETATM, etc)

_ characterize the interface (sugar/phosphate/base - protein distances, water contacts, etc)

_ use the 3DNA program [1] for NA structure description that gives exhaustive data in easily parsable Json format

_ rearrange the data per nucleotide (eg. “nucl 5 to 15 make a stem-loop” → “nucl 5 is at position 1 in an 11-nucl stem-loop”)

_ optionaly: cut the NA structures into fragments (trinucleotides in our example)

[1] X-J Lu & WK Olson. 3DNA: a software package for the analysis, rebuilding and visualization of three‐dimensional nucleic acid structures. Nucleic Acids Research (2003) 31(17), 5108-21

The output is a description, at both the structural and (tri)nucleotide levels, of the full set of protein-bound NAs from the PDB, in a single Json file.

--------------------------------------------------------------------------
Requests to select structures / fragments
--------------------------------------------------------------------------
[to be written in march 2019]

--------------------------------------------------------------------------
Requests to compute statistics on the database
--------------------------------------------------------------------------
[to be written in march 2019]
