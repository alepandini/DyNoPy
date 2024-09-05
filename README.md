# DyNoPy (Dynamics based Network cOmparisons in Python)

Welcome to DyNoPy!!!!

DyNoPy, for your favourite protein, will allow you to combine coevolution matrix with data from Molecular Dynamics (MD) simulations to find important residues that might be contributing to function.

So before you start you need the following:

Input files:
1. PDB file used for MD simulation, residue numbering, atom count should match the trajectory file.
2. Trajectory file, no. of atoms shoul match the PDB file.
3. FASTA sequence, should match the length and sequence in the pdb file

... and most importantly DyNoPy installed on your workhorse server/PC along with the dependency softwares listed below.

## **What does it do?**
This package will help you generate:
1. Coevolution matrices from FASTA sequences
2. Pairwise non-bonded interaction energies from MD trajectories
3. Pairwise residue contributions to functional motions
4. Pairwise combined covevolution and dynamical score

## **Data analysis and visualisation**
Use the sample Jupyter notebooks in the jupyter/ folder for data visualisation and analysis

---

## How to Install

Below is the list of softwares and the database required to run DyNoPy. install.sh script will do the work for you, but if you are interested in the details, you can find them below:

Please download the latest version of uniclust database:
Link to [uniprot databases](http://gwdu111.gwdg.de/~compbiol/uniclust/)

Sequence alignment
- [hhsuite](https://github.com/soedinglab/hh-suite)
Coevolution matrix (any would be fine, we recommend CCMpred for its speed and accuracy)
- [CCMpred](https://github.com/soedinglab/CCMpred)
- [Metapsicov](https://github.com/psipred/metapsicov)
- [plmDCA](https://github.com/pagnani/PlmDCA)
- [ambertools](https://ambermd.org/AmberTools.php)
---

#### **Python requirements**
- Python3.6, numpy 1.15, matplotlib
- R & igraph
	
### **Installation instructions**

Download the code to your favourite directory
```
cd /home/username/myfavdir
git clone git@github.com:alepandini/DyNoPy.git
```

This is to make your life easy. Add the following lines to your:
```
vim .bashrc
```

For ubuntu to figure out where DyNoPy executables are
```
export DYNOPY="/home/username/myfavdir/DyNoPy"
export PATH=$PATH:$DYNOPY/bin
```
for python to figure out where DyNoPy python code is
```
export PYTHONPATH=$PATH:$DYNOPY/bin
```
... in your favourite dir i.e. `/home/username/myfavdir` run the following command to install cmake, hh-suite, and CCMpred
```
dynopy_installer_00.sh
```
.bashrc settings for dependcy softwares. If these variables are not set, `dyno_coevolution.py` will complain alot

HHBLITS
```
export HHLIB="/home/username/myfavdir/hh-suite"
export PATH="$PATH:$HHLIB/build/bin"
```
CCMPred
```
export CCMPRED_HOME="/home/username/myfavdir/CCMpred"
export PATH="$PATH:$CCMPRED_HOME/build/bin"
export PATH="$PATH:$CCMPRED_HOME/scripts"
```
After you download and install ambertools 
```
export AMBERHOME="/home/username/myfavdir/amber/amber18"
export PATH="$PATH:$AMBERHOME/bin"
```

Now that you have installed DyNoPy, head to [tutorials](https://github.com/alepandini/DyNoPy/wiki) on how to use DyNoPy.
