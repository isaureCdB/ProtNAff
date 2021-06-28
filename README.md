# ProtNAff
Create filters to select Protein - Nucleic acids structures from PDB files, and build libraries of protein-boud RNA fragments

--------------------------------------------------------------------------
ProtNAff is a pipeline to:
--------------------------------------------------------------------------
1. Clean-up and parse all NA-protein structures from the PDB into ensembles of small information units in a single file

2. Search for sets of NA-protein structures with highly cutomisable combinations of criteria

3. Create RNA/DNA 3D fragment libraries extracted from those sets of structures.

4. Perform statistics on customised features of such libraries

1 is the necessary first step, 2 and 3 can be done independantly, 4 must be done after 3.
The output of 1 for the PDB at a certain time will soon be downloadable from the LORIA website, allowing to do directly 2 and/or 3(+/-4) directly.

--------------------------------------------------------------------------
### Installation
--------------------------------------------------------------------------

For the installation come [here](./INSTALLATION.md)

--------------------------------------------------------------------------
### Testing and Examples
--------------------------------------------------------------------------

There are several notebooks to help you to understand ProtNAff:

- The [example notebook](./example/example.ipynb) is a notebook that will
help you to create a small database and your first fragment library. Add the
end of the notebook you will do a small graph, to know is the installation is
correct you can compare this graph to the one in the next notebook.

- The [test notebook](./example/test.ipynb) is a small notebook that create
the graph that you need to obtain at the end of the previous notebook, if both
graph are the same, then everything is fine.

- The [data_protnaff notebook](./data_protnaff.ipynb) is a notebook that will
create the JSONs containing the data used for the paper.

- The [figures_protnaff notebook](./figures_protnaff.ipynb) is the notebook that
uses the JSONs from `data_protnaff.ipynb` to create the figures from the paper.

--------------------------------------------------------------------------
### Work in progress
--------------------------------------------------------------------------
The main idea of ProtNAff is to provide a highly versatile pipeline to cover as many usages as possible.
This is intended to be a dynamic collaborative work:

_ If you need a specific feature that you can't find or don't know how to add in the current pipeline, please contact us, and we will do our best to include it!

_ If you added some feature that you think can be useful to others, please fill free to propose a pull request.
