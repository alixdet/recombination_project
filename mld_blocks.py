"""
Long project
SMILE unit

Ref : 
- Kerdoncuff, Lambert & Achaz, 2020
"""

import os
import random

from collections import Counter


"""Generate sequences
should be done with msprime but doesn't work
"""
len_seq = 500
alignment = []

random.seed(1234)

alignment = [''.join(random.choices(('A', 'C', 'G', 'T'), k = len_seq)) \
                for i in range(4)]

# same lengths by construction but otherwise should be verified
for seq in alignment[1:]:
    assert len(seq) == len(alignment[0]), 'Sequences dont have same length'


"""Alignment
to cleanly import an existing alignment,
see https://www.tutorialspoint.com/biopython/biopython_sequence_alignments.htm

or, to align using MUSCLE codewrapper:
https://biopython.org/docs/1.75/api/Bio.Align.Applications.html

DOESNT WORK
"""
#from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline

muscle_exe = r'/home/alix/.softwares/muscle_v5.0.1428_linux'
in_file = r'unaligned.fasta'
out_file = r'aligned.fasta'

muscle_cline = MuscleCommandline(muscle_exe,
                                 input=in_file,
                                 out=out_file)

#print(muscle_cline)
assert os.path.isfile(muscle_exe), 'MUSCLE 5 executable missing'
# #Non-zero return code 1 from ... , message 'Invalid command line'
# stdout, stderr = muscle_cline()


"""Get list of intervals
"""

def sign_polymorphisms(bases :list):
    """
    4-gametes test function

    From a list of 4 elements, returns positions of pairwise polymorphisms,
    ("significant"), else False
    """
    assert len(bases) == 4, "list of 4 elments required"

    counter = Counter(bases)
    pos = []

    if list(counter.values()) != [2, 2]:
        return False

    for base in counter.keys():
        pos.append([index for index, value in enumerate(bases) \
            if value == base])
    return pos


########################################
##### First way ########################
########################################

# first 150 polymorphic sites
intervals = []
polymorphism_met = 0
i = 0

first_polymorphism = False  # store position of first polym

# should catch IndexOutOfRange error
while polymorphism_met < 10: # P = 150 in paper
    polym = sign_polymorphisms([seq[i] for seq in alignment])
    
    if polym:
        # should consider all polymorphisms ? i.e. not just pairwise ?
        polymorphism_met += 1
        if not first_polymorphism: # i.e. no polym previously met
            left_border = i
            first_polymorphism = polym
        else :
            if polym != first_polymorphism: # incompatibility found !
                intervals.append([left_border, i])
                # a break could start chopping right here !
    i += 1


########################################
### Second way (BETTER) ################
########################################
intervals = []

# first get alignment length required to meet given number of polymorphisms
P = 10
polymorphism_met = 0

i = 0
while polymorphism_met < P:
    # should consider all polymorphisms ? i.e. not just pairwise ?
    if sign_polymorphisms([seq[i] for seq in alignment]):
        polymorphism_met += 1
    i += 1
requi_len = i


for i in range(requi_len):
    first_polym = sign_polymorphisms([seq[i] for seq in alignment])

    if first_polym : # browses genome from first polym
        for j in range(i + 1, requi_len):
            second_polym = sign_polymorphisms([seq[j] for seq in alignment])
            if second_polym and second_polym != first_polym: # ugly way
                intervals.append((i,j))


"""Shortening list of intervals

intervals list is sorted and defined the way it is above
-> does not work otherwise

try with yield statements
"""

intervals_rev = intervals[::-1]

chopped_intervals = [intervals_rev[0]]

for i in range(1, len(intervals_rev)):
    if intervals_rev[i][1] > intervals_rev[0][0]:
        chopped_intervals = [(intervals_rev[0][0], intervals_rev[1][1])]


    
"""MLD breakpoints
"""

