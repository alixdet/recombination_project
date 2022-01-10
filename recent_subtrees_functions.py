"""
Long project
SMILE unit

Ref : 
- Kerdoncuff, Lambert & Achaz, 2020
"""

from alignment_to_mld_functions import *

def avg_mutation_rate(alignment :list):
    """
    Counts number of POLYMORPHIC sites on an alignment.
    -> no distinction on whether mutation carried by 1 or several sequences
    Used to estimate the age of the most recent common ancester.

        Parameters:
            alignment (list): !should be a compatible segment
        Outputs:
            float: #polymorphisms / len(aln) x #seq
    """
    count = 0
    for bases in alignment:
        if len(set(bases)) > 1:
            count += 1

    return count / (len(alignment) * len(alignment[0]))


def is_recent(alignment :list, th :float):
    """
    Whether the TMRCA of a coalescent tree is likely recent.
    Calls avg_mutation_rate().

        Parameters:
            alignment (list): !should be a compatible segment
            th (float): threshold !!! TO BE DETERMINE !!!
        
        Output:
            boolean
    """
    # PROBABLY USELESS
    pass


