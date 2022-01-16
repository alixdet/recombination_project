# Recombination project

Master 2 Bio-informatics, Université de Paris & *Stochastic Modelisation for the Inference of Life Evolution unit*, Collège de France


Reference :
Testing for population decline using maximal linkage disequilibrium blocks, Kerdoncuff, Lambert & Achaz, *Theoretical Population Biology*, 2020

# Motivation
This project aims to continue the work conducted by the lab on population dynamics and on the use of *Maximum Linkage Disequibrilium (MLD)* blocks to infer recent variations in a population. The goal is to detect subtrees within a recombination-free block, whom time of the most recent common ancestor is below is thershold, i.e. recent subtrees.

# Usage
## Installation
Consider installing the requirements with *conda* from the main directory:
```
$ cd recombination_projet
$ conda env create -f environment.yml
$ conda activate recombination_projet
```

## Run
The program requires as input an sequence alignment in a single fasta file.

It can then be run with the command in the terminal:
```
$ python3 main_aln_to_polym-ratio.py -a aln/aln_100k.fa -s 0 -e 1000
```
It takes as argument the path to the file (-a), the start (-s) and the end (-e) positions on the alignment.
The output is the distribution of the ratio of polymorphic sites within MLD blocks.







