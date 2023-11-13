VAR_TMP=`python -c "import os, sys; print(os.path.abspath(os.path.dirname(sys.argv[1])))" $0`
echo "export PROTNAFF=$VAR_TMP" >> ~/.bashrc
echo 'export PYTHONPATH=$PROTNAFF:$PYTHONPATH' >> ~/.bashrc
