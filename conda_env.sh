#! /bin/sh
sudo apt install libcairo2-dev 
conda update -n base -c defaults conda
conda deactivate
conda remove --name dynopy --all
conda env create -f environment.yml
conda activate dynopy
