#! /bin/sh
conda deactivate
conda remove --name dynopy --all
conda env create -f environment.yml
conda activate dynopy
