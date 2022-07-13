# mjolnir

Pipeline tools

Link to uniprot databases
	- http://gwdu111.gwdg.de/~compbiol/uniclust/2022_02/

Sequence alignment
- hhsuite https://github.com/soedinglab/hh-suite

Coevolution matrix
- ccmpred https://github.com/soedinglab/CCMpred
- metapsicov https://github.com/psipred/metapsicov
- plmDCA https://github.com/pagnani/PlmDCA

Secondary structure calculation

PsiPRED:
- https://github.com/psipred/psipred

# DyNoPy (Dynamics based Network cOmparisons in Python)


## **What does it do?**
This package will help you generate:
1. Coevolution matrices from FASTA sequences
2. Pairwise non-bonded interaction energies from MD trajectories 
3. Pairwise residue contributions to functional motions
4. Pairwise combined covevolution and dynamical score

#### **Requirements**
	- Python3.6, numpy 1.15, matplotlib
	- R & igraph
        - ambertools18
        - hh-suite
	- CCMpred
	
## **Installation instructions**

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

## **Developer area**

### **Progress**
	- Convert tools to wrappers for class calls
	- Clean up the code
- **Co-evolution analysis**
    - *Features to add*
       - [ ] check for ccmpred/hhblits libraries
- **RIE**
    - *Features to add*
       - [ ] check for ambertools
- **Dependency files**
    - [ ] mdp files
    - [ ] R scripts
### GitHub cheat sheet
Some tips and tricks of GitHub

##### Configure git 
```
git config --global user.name "Name"
git config --global user.email "email id"
git config --global core.editor "vim"
```
#### Setting up the project 
```
mkdir DyNoPy
cd DyNoPy/
git init
git remote add DyNoPy git@github.com:alepandini/DyNoPy.git
#check if it is setup properly
git remote -v
git pull DyNoPy master
```

##### add+commit+push
```
cd DyNoPy && git add . && git commit && git push DyNoPy master
```
