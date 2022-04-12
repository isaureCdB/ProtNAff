# ProtNAff
Create filters to select Protein - Nucleic acids structures from PDB files, and build libraries of protein-bound RNA fragments

--------------------------------------------------------------------------
ProtNAff is a pipeline to:
--------------------------------------------------------------------------
1. Clean-up and parse all NA-protein structures from the PDB into ensembles of small information units in a single file.

2. Search for sets of NA-protein structures with highly cutomisable combinations of criteria.

3. Create RNA/DNA 3D fragment libraries extracted from those sets of structures.

4. Perform statistics on customised features of such libraries.

Step 1 is the necessary first step, steps 2 and 3 can be done independantly, step 4 can be done only after step 3.
The output of 1 for the PDB at a certain time will soon be downloadable from the LORIA website, allowing to do steps 2, 3 and/or 4 directly.

--------------------------------------------------------------------------
### Installation
--------------------------------------------------------------------------

Installation instructions are [here](./INSTALLATION.md)

--------------------------------------------------------------------------
### Testing and Examples
--------------------------------------------------------------------------

There are several notebooks to help you to understand ProtNAff:

- The [example notebook](./example/example.ipynb) notebook helps 
you create a small database and your first fragment library. At the
end of the notebook, a graph is creates to check if the installation is
correct, by comparing this graph to the one in the next notebook.

- The [test notebook](./example/test.ipynb) notebook creates the same graph 
as should be obtained at the end of the example notebook: if both
graphs are identical, the installation went fine.

- The [data_protnaff notebook](./data_protnaff.ipynb) notebook creates 
the JSON files containing the data uses for analyses in the protNAff paper.

- The [figures_protnaff notebook](./figures_protnaff.ipynb) notebook uses 
the JSON files created by `data_protnaff.ipynb` to create the figures in the paper.

- The [figures_dna_protnaff notebook](./figures_dna_protnaff.fr) notebook creates
the same figures as previously but for DNA instead of RNA.

--------------------------------------------------------------------------
### Work in progress
--------------------------------------------------------------------------
The main idea of ProtNAff is to provide a highly versatile pipeline to cover as many usages as possible.
This is intended to be a dynamic collaborative work:

_ If you need a specific feature that you can't find or don't know how to add in the current pipeline, please contact us, and we will do our best to include it.

_ If you added some feature that you think can be useful to others, please fill free to propose a pull request.
