"""
Long project
SMILE unit

Ref : 
- Kerdoncuff, Lambert & Achaz, 2020
"""

import argparse
import os
import sys
from os import listdir, getcwd
from os.path import isfile, join

from alignment_to_mld_functions import *
from recent_subtrees_functions import *


def get_arguments():
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                        "{0} -h"
                                        .format(sys.argv[0]))
    parser.add_argument('-d',
                        dest='dir',
                        required=True,
                        help="Directory with aligned sequences")
    return parser.parse_args()


def get_bases_files(bfs :list, pos :int, offsets :list):
    """
    !: filters on 'N's and case

        Parameters:
            bfs (list): buffer files
            pos (int): position on the alignement
            offests (list of int): offsets on different alignement files
                to decount for the reading
        
        Output:
            list: bases on the position
    """
    assert len(bfs) == len(offsets), "Buffers or offset values missing."

    bases = []
    for i in range(len(bfs)):
        bfs[i].seek(offsets[i] + pos + pos//60)
        bases.append(bfs[i].read(1).upper())

    if 'N' not in bases:
        return bases


def main():
    args = get_arguments()

    wd = getcwd() + args.dir
    seqfiles = [os.path.join(wd, f) for f in os.listdir(wd) if \
        os.path.isfile(os.path.join(wd, f))]
    opened = []

    bfs = ['bfs'+ str(i) for i, f in enumerate(seqfiles)]

    # Cannot be done within a function as the files need to remain open
    # during the queries

    # used to decount first line
    offsets = ['NaN'] * len(seqfiles)

    for i, seqfile in enumerate(seqfiles):
        bfs[i] = open(seqfile, "r")
        opened.append(bfs[i])
        bfs[i].readline()
        offsets[i] = bfs[i].tell()
    
    nb_seq = len(opened)

    l = 2000

    ###
    # Cutting in MLD
    ###
    # get list of incompatible polymorphism along the genome
    ip = incompatible_polym([get_bases_files(bfs, i, offsets) \
        for i in range(l)])

    print(ip)
    print('---')

    chopped = chop_list(ip)
    mld = place_mld_breakpoints(chopped)
    mld.sort()
    #print(mld[:10])

    ###
    # Date the trees
    ###
    mutation_rates = []

    # EDGE EFFECT ! : last block not taken into account
    # -> add total length to list
    # TODO: find out how to grap total length without computing entire file
    mld = [0] + mld + [l]

    for i in range(len(mld) - 1):
        count = 0
        for pos in range(mld[i],mld[i+1]):
            segment_len = mld[i+1] - mld[i]
            bases = get_bases_files(bfs, pos, offsets)
            if len(set(bases)) > 1:
                count += 1
        mutation_rates.append(count / (segment_len * nb_seq))

    print(len(mld))
    print(mld)
    print('---')
    print(len(mutation_rates))
    print(mutation_rates)


    [f.close() for f in opened]


if __name__ == "__main__":
    main()

# to try 
# python3 main_alignment_to_mld.py -d /Hironda_rustica_4/