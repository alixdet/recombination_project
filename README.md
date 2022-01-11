# Recombination project

Master 2 Bio-informatics, Université de Paris & SMILE unit, Collège de France

1 month project under the supervision of Thomas Forest, Guillaume Achaz & Amaury Lambert.

Reference :
Testing for population decline using maximal linkage disequilibrium blocks, Kerdoncuff, Lambert & Achaz, *Theoretical Population Biology*, 2020

# Motivation
This project aims to continue the work conducted by the lab on population dynamics and on the use of *Maximum Linkage Disequibrilium (MLD)* blocks to infer recent variations in a population. The goal is to detect subtrees within a recombination-free block, whom time of the most recent common ancestor is below is thershold, i.e. recent subtrees.

# Usage
## Installation
Consider installing the requirements with *conda* from the main directory:
```
$ cd mld_recent_subtrees
$ conda env create -f environment.yml
$ conda activate mld_recent_subtrees
```

## Run
The program requires as input a directory containing sequences aligned in fasta files. The files need to be in standard fasta format, the first line starting with '>' followed by the sequenced ID, and the sequences on 60 characters-long lines.

It can then be run with the command in the terminal:
```
$ python3 main_alignment_to_mld.py -d /dir_with_files/
```







