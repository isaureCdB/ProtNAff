# Installation of ProtNAff

To install ProtNAff you need to be on Linux, it will not work on Windows and on MacOS.

1. First you need to clone the repository using Git. To check if you have Git on your computer you can execute this command in a terminal:

`git --version`

If you have a git version you can continue, else you need to install Git (https://git-scm.com/downloads)

2. To clone the repository you can use this command (first go to the directory where you want to have ProtNAff):

`git clone https://github.com/isaureCdB/ProtNAff.git`

This will create a ProtNAff folder.

3. You will also need miniconda3 (https://docs.conda.io/en/latest/miniconda.html), as you will use several environement to switch from python2 to python3.

When miniconda is installed you will be able to use both environement of ProtNAff, to do so you will need to lauch those lines into the ProtNAff folder:

`conda env create -f protnaff_environment.yml`

`conda env create -f attract_environment.yml`

4. You need 1 more tool which is jq, on Ubuntu you can run the following command:

`sudo apt install jq`

4. To finish you need to add several variable of environement. To do so you have add them in your .bashrc.
We created a bash file that will do it for you, run the following command:

`bash add_env_var.sh`

5. Then you need to source your .bashrc or open a new terminal.

6. The how_to_use.ipynb file is a jupyter notebook file, it will help you to understand how to use the tool.
To run it you need to run 2 command lines into the example folder:

conda activate protnaff

jupyter notebook

When this is done a web page should open, it is the notebook, to execute a cell you have to press Ctrl+Enter.

If you have any question, do not hesitate to contact us !
