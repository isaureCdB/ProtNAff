# User manual for ProtNAff

Protein-bound Nucleic Acid filters and fragment libraries is a tool designed for bioinformatic.
It is made to create filters to select Protein - Nucleic acids structures from PDB files, and build libraries of protein-bound RNA fragments.

This document will provide explanations of how to use ProtNAff and what you should be able to do with it.

--------------------------------------------------------------------------
To sumarize a bit, ProtNAff is a pipeline to:
--------------------------------------------------------------------------
1. Clean-up and parse all NA-protein structures from the PDB into ensembles of small information units in a single file.

2. Search for sets of NA-protein structures with highly customisable combinations of criteria.

3. Create RNA/DNA 3D fragment libraries extracted from those sets of structures.

4. Perform statistics on customised features of such libraries.

Step 1 is the necessary first step, steps 2 and 3 can be done independantly, step 4 can be done only after step 3.
The output of 1 for the PDB at a certain time will soon be downloadable from the LORIA website, allowing to do steps 2, 3 and/or 4 directly.

### Creation of the structures database

The first step is to create the file named `structures.json`. This file is a JSON containing information extract from x3DNA, but also some information computed by ProtNAff. The example usages are given lower.

### Usage of filters

The filters are the most customizable part, indeed it is up to you to create the filters you need.
Example of filters are given in the paper, but there are also other example in the `filters` folder.

Some of them are python scripts as `filter_no_modified.py`, which is returning a JSON file, or
`filter_ss.py`, which is printing the single-stranded nucleotides per PDB id.

The detailed information contained in the `structures.json` are:

* Per PDB structure:
  - Experimental method;
  - Resolution (if X-ray);
  - Name of NA chains and protein chains;
  - Number of models;
  - Name of cofactors.

* Per nucleic acid chain:
  - Position of breaks within the backbone;
  - Sequence.

* Per nucleotide:
  - All H-bonds made with the protein, with for each H-bond (i) the amino-acid type and atom, (ii) the nucleotide sub-part, (iii) the H-bond distance, and (iv) the fact that 3DNA would consider the H-bond donor or acceptor as "questionable"
  - All H-bonds made with another nucleic acid, with (i) the position of the other nucleotide in the sequence (n-2, n-1, n+1, n+2 or other), (ii) the nucleotide sub-part, (iii) the H-bond distance, and (iv) the fact that 3DNA would consider the H-bond donor or acceptor as "questionable" according to DSSR;
  - Total number of H-bonds with protein for each sub-part, with a 0.5 weighting when questionable;
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


All those information can be used to filter the structures or the fragment depending of what you want.
It is also possible tu use the 3DNA output directly, as we are doing in the `filter_hairpin.py`.

### Creation of fragment libraries

Using the sets of structure made in the previous step, it is possible to do fragment libraries.
The clustering can be done using 2 methods. The first one is the method named "fastclust", and the second one which is more accurate but slower is the radius method.

### Statistics

The statistics are completely up to the user. Nothing is provided excepted the examples from the paper in the notebook named below.

--------------------------------------------------------------------------
### Installation
--------------------------------------------------------------------------

Installation instructions are [here](./INSTALLATION.md)

**In addition, you can run the ProtNAff filtering and clustering web server by clicking [here](https://colab.research.google.com/github/isaureCdB/ProtNAff/blob/master/filtering-clustering.ipynb).**

--------------------------------------------------------------------------
### Testing and Examples
--------------------------------------------------------------------------

There are several notebooks to help you to understand ProtNAff:

- The [example notebook](./example/example.ipynb) helps
you create a small database and your first fragment library. At the
end of the notebook, a graph is creates to check if the installation is
correct, by comparing this graph to the one in the next notebook.

- The [test notebook](./example/test.ipynb) creates the same graph
as should be obtained at the end of the example notebook: if both
graphs are identical, the installation went fine.

- The [data_protnaff notebook](./data_protnaff.ipynb) creates
the JSON files containing the data uses for analyses in the protNAff paper.

- The [figures_protnaff notebook](./figures_protnaff.ipynb) uses
the JSON files created by `data_protnaff.ipynb` to create the figures in the paper.

- The [figures_dna_protnaff notebook](./figures_dna_protnaff.ipynb)creates
the same figures as previously but for DNA instead of RNA.


- The [filtering and clustering notebook](./filtering-clustering.ipynb) shows you how to perform custom filtering and clustering.
In the notebook, you can select and run a custom filter and clustering method among the provided examples. You can also write your own filter or clustering method.
**You can also run this notebook as a web server in Google Colab by clicking [here](https://colab.research.google.com/github/isaureCdB/ProtNAff/blob/master/filtering-clustering.ipynb). This does not require ProtNAff to be installed.**


### Description of the main scripts

Here the main scripts are decribed:

* `create_database.sh` the script to create the `structures.json`
  - input: a list of pdb ids, and rna/dna depending of what nucleic acids you are working on
  - ouput: the `structures.json` file

* `create_frag_library.sh` the script to create the fragment libraries
  - input: rna/dna depending of what nucleic acids you are working on, the `structures.json` is needed
  - output: `fragments.json` and `fragments_clust.json`

Some useful scripts:

* `create_frag_library/npy2pdb.py` this script is converting the npy matrix into pdb files.
  - input: the npy matrix, and the pdb template, corresponding to the atoms
  - output: a string output containing all the lines corresponding to the pdb information, it can be cast into a pdb file

* `create_frag_library/pdb2npy.sh` this script convert the pdb files into a npy matrix
  - input: a list of pdb file to convert to npy, --list, --outp the name of the matrix you will obtain
  - output: a npy matrix

* `create_frag_library/reduce.py` this script is reducing to coarse grained the inputs
  - input:
  - output: 

--------------------------------------------------------------------------
### Work in progress
--------------------------------------------------------------------------
The main idea of ProtNAff is to provide a highly versatile pipeline to cover as many usages as possible.
This is intended to be a dynamic collaborative work:

_ If you need a specific feature that you can't find or don't know how to add in the current pipeline, please contact us, and we will do our best to include it.

_ If you added some feature that you think can be useful to others, please fill free to propose a pull request.
