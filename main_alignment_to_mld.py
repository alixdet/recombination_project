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

from four_gametes import *


def read_fasta(seq_path :str):
    seq_list = []
    with open(seq_path, 'r') as seq_w:
        for line in seq_w:
            if not line.startswith('>'):
                seq_list.append(line.rstrip('\n'))
    return seq_list


def get_arguments():
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                        "{0} -h"
                                        .format(sys.argv[0]))
    parser.add_argument('-d',
                        dest='dir',
                        required=True,
                        help="Directory with aligned sequences")
    # parser.add_argument('-o',
    #                     dest='output_file',
    #                     type=str,
    #                     default=os.curdir + os.sep + "contigs.fasta",
    #                     help="Output contigs in fasta file")
    return parser.parse_args()


def main():
    args = get_arguments()

    wd = getcwd() + args.dir
    seqfiles = [os.path.join(wd, f) for f in os.listdir(wd) if \
        os.path.isfile(os.path.join(wd, f))]

    for seqfile in seqfiles:
        with open(seq, 'r') as seq_file:
            for read in read_fasta(args.reads_file):
                d = match(read, N, fm_index)
                of.write(str(pos(d[0], d[1], pos_bwt)))

if __name__ == "__main__":
    main()
