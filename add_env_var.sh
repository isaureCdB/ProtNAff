echo 'export MY_CONDA=~/miniconda3' >> ~/.bashrc
VAR_TMP=`python -c "import os; print(os.getcwd())"`
echo "export PROTNAFF=$VAR_TMP" >> ~/.bashrc
echo 'export PYTHONPATH=$PROTNAFF:$PYTHONPATH' >> ~/.bashrc
