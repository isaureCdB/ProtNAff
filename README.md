# User manual for ProtNAff

Protein-bound Nucleic Acid filters and fragment libraries (ProtNAff) is a tool to create filters selecting structures of Protein - Nucleic acids complexes from the PDB and to build libraries of protein-bound RNA fragments.

This document explains how to use ProtNAff and what can be done with it.

--------------------------------------------------------------------------
To sumarize, ProtNAff is a pipeline to:
--------------------------------------------------------------------------
1. Clean-up and parse all NA-protein structures from the PDB into ensembles of small information units in a single file.

2. Search for sets of NA-protein structures with highly customisable combinations of criteria.

3. Create RNA/DNA 3D fragment libraries extracted from those sets of structures.

4. Perform statistics on customised features of such libraries.

Step 1 is necessary, steps 2 and 3 can be done independantly, step 4 can be done only after step 3.
The output of step 1 (for the PDB at a certain time) can now be downloaded from https://zenodo.org/record/6475637#.YmK3MFxByV4 , allowing to do steps 2, 3 and/or 4 directly.

### Creation of the structures database (step 1)

The first step is to create `structures.json`, a JSON containing informations per 3D structure, either extracted from x3DNA outputs or computed by ProtNAff directly. Examples of usage are given below.

### Usage of filters (step 2)

The filters are the most customizable part: you can create the filters you need.
Examples of filter are given in the ProtNAff paper, and more examples are given in the `filters` folder. Especially, look at the
explanation_filters.ipynb notebook for a detailed explanation on how to build filters.

Some filters are python scripts, such as `filter_no_modified.py`, which returns a JSON file, or
`filter_ss.py`, which prints the single-stranded nucleotides per PDB id. Others are jupyter-notebooks, such as filtering-clustering.ipynb. The earliers must be runon your machine after installation of the proper conda environement, while the laters can be run through a collab session (see *INSTALLATION* section below).

The detailed information contained in the `structures.json` are:

* Per PDB structure:
  - Experimental method;
  - Resolution (if X-ray);
  - Name of NA chains and protein chains;
  - Number of models (for NMR);
  - Name of cofactors.

* Per nucleic acid chain:
  - Position of breaks within the backbone;
  - Sequence.

* Per nucleotide:
  - All H-bonds made with the protein, with for each H-bond (i) the amino-acid type and atom, (ii) the nucleotide sub-part (phosphate group/sugar/base), (iii) the H-bond distance, and (iv) the fact that 3DNA would consider the H-bond donor or acceptor as "questionable"
  - All H-bonds made with another nucleic acid, with (i) the position of the other nucleotide in the sequence (n-2, n-1, n+1, n+2 or other), (ii) the nucleotide sub-part, (iii) the H-bond distance, and (iv) the fact that 3DNA would consider the H-bond donor or acceptor as "questionable" according to DSSR;
  - Total number of H-bonds with protein for each sub-part (phosphate group/sugar/base), with a 0.5 weighting for questionable H-bonds according to 3DNA;
  - Base-pairing types it is involved in;
  - Initial name of the residue in the PDB file (if canonized residue);
  - Minimal distance of each sub-part to the protein and to cofactors, if < 5 \AA;
  - Parts that had missing atoms in the initial PDB file;
  - Secondary structure (terminal single-stranded parts, hairpin loop, internal loop, junction, double-stranded);
  - Presence of a stacking interaction with nucleotides at position n-2, n-1, n+1, n+2 or any position in sequence.

* Per fragment:
  - Name of the PDB structure it is extracted from;
  - Model index (for multi-model PDB structures);
  - Name of the PDB chain;
  - Residue indices in the PDB file;
  - Initial sequence;
  - In which part of which nucleotide were atoms missing (if any);
  - If the fragment is a cluster prototype for the different clustering thresholds (0.2\AA, 1\AA ~and 3\AA ~in the current implementation);
  - Index of the cluster it belongs to, for the different thresholds.

All those information can be used and combined to filter the structures or the fragments and create a set suitable for a given application.
It is also possible tu use the 3DNA outputs directly, as we are doing in the `filter_hairpin.py`.

### Creation of fragment libraries (step 3)

This step creates a 3D fragment library from the set of structures created in the previous step.
The clustering creates clusters of fragments that are at a maximal RMSD from the cluster center. This can be done by two methods: 
- "fastclust", which is fast but non deterministic (dependant on the fragments order) and does not minimize the total number of clusters
- "radius", which is slow but deterministic and minimize the number of clusters.

### Statistics (step 4)

ProtNAff allows to run all kinds of statistics on the structures database and on the fragment libraries. Examples of statistics from the protNAff paper are provided in the notebooks named below (*Testing and Examples* section).

--------------------------------------------------------------------------
### Installation
--------------------------------------------------------------------------

Installation instructions are [here](./INSTALLATION.md)

**Alternatively, you can run the ProtNAff filtering and clustering web server by clicking [here](https://colab.research.google.com/github/isaureCdB/ProtNAff/blob/master/filtering-clustering.ipynb).**

--------------------------------------------------------------------------
### Testing and Examples
--------------------------------------------------------------------------

There are several notebooks to help you to understand ProtNAff:

- The [example notebook](./example/example.ipynb) helps
you create a small database and your first fragment library. At the
end of the notebook, a graph is creates to check if your installation is
correct, by comparing this graph to the one in the next notebook.

- The [test notebook](./example/test.ipynb) creates the same graph
as should be obtained at the end of the example notebook: if both
graphs are identical, the installation went fine.

- The [data_protnaff notebook](./data_protnaff.ipynb) creates
the JSON files containing the data used for analyses in the protNAff paper.

- The [figures_protnaff notebook](./figures_protnaff.ipynb) uses
the JSON files created by `data_protnaff.ipynb` to create the figures in the protNAff paper.

- The [figures_dna_protnaff notebook](./figures_dna_protnaff.ipynb) creates
the same figures as previously but for DNA instead of RNA.

- The [filtering and clustering notebook](./filtering-clustering.ipynb) shows you how to perform custom filtering and clustering.
In the notebook, you can select and run custom filter and clustering methods among the provided examples. You can also write your own filter or clustering method.
**You can also run this notebook as a web server in Google Colab by clicking [here](https://colab.research.google.com/github/isaureCdB/ProtNAff/blob/master/filtering-clustering.ipynb). This does not require ProtNAff to be installed.**


### Description of the main scripts

The main protNAff scripts:

* `create_database.sh` to create the `structures.json` file.
  - input: a user-given list of pdb ids of NA-protein complexes and the type of nucleic acids you work on (rna/dna). Optional: a threshold for NA-protein/cofactor contact distance'.
  - ouput: the `structures.json` file
  - usage: ./create_database.sh <list of PDBs> --rna/dna -t <float>

* `create_frag_library.sh` the script to create the fragment libraries
  - input: `structures.json`, and the type of nucleic acids you work on (rna/dna)
  - output: `fragments.json` and `fragments_clust.json`
  - usage: ./create_frag_library.sh --rna/dna

Some useful scripts:

* `create_frag_library/npy2pdb.py` converts an npy matrix of 3D coordinates into pdb files.
  - input: the npy matrix, and a pdb template (with the same atoms in same order). Optional: an index or a list of indices of the structures to be converted into PDB.
  - output: a multi-model pdb file
  - usage: ./create_frag_library/npy2pdb.py filename.npy template.pdb [--list <list of indices>] [--index <integer>] > filename.pdb

* `create_frag_library/pdb2npy.py` converts pdb files into a npy matrix.
  - input: a text file containing a list paths to the pdb files to be converted into npy 
  - output: a npy matrix
  - usage: create_frag_library/pdb2npy.py filename.list --list --outp filename.npy

* `create_frag_library/reduce.py` this script reduces all-atom rna/dba pdb into ATTRACT coarse-grained representation

--------------------------------------------------------------------------
### Work in progress
--------------------------------------------------------------------------
The main idea of ProtNAff is to provide a highly versatile pipeline to cover as many usages as possible.
This is intended to be a dynamic collaborative work:

_ If you need a specific feature that you can't find or don't know how to add in the current pipeline, please contact us, and we will do our best to include it.

_ If you added some feature that you think can be useful to others, please feel free to propose a pull request.
