# DyNoPy
Dynamics based Network cOmparisons in Python
#### **What does it do?**
This package will help you extract
	- Coevolution matrix from a FASTA sequence
	- Extract pairwise non-bonded interaction energies from MD trajectories 
	- Identify residue pairs significant for functional motion
	- Identify communities of residues in a network and rank them 
	- Compare networks (how to treat deletions??)

#### **Requirements**:
	- Python3.6, numpy 1.15, matplotlib
	- R & igraph 
	

#### **Progress**:
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
#### git cheat sheet
##### Configure git 
```
git config --global user.name "Name"
git config --global user.email "email id"
git config --global core.editor "vim"

mkdir DyNoPy
cd DyNoPy/
git init
git remote add DyNoPy git@github.com:alepandini/DyNoPy.git
#check if it is setup properly
git remote -v
git pull DyNoPy master
```
##### add+commit+git
```
cd DyNoPy && git add . && git commit && git push DyNoPy master
```
