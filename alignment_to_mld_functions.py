"""
Long project
SMILE unit

Ref : 
- Kerdoncuff, Lambert & Achaz, 2020
"""

from collections import Counter


def info_polymorphism(bases :list):
    """
    Returns whether the polymorphism is informative,
    i.e. 2 alleles at least twice.

    Parameters: 
            bases (list): 4 or more bases

    Returns:
            boolean: polymorphism is informative
    
    TODO: consider case of random mutation on informative polymorphism
        when # seq >= 5
    """
    assert len(bases) >= 4, "4 or more elments required"

    # maybe find a way w/o using Counter
    values = Counter(bases).values()

    return len(values) >= 2 and 1 not in values

def polym_positions(alignment :list):
    """
    Looks for informative polymorphism on an alignment, # seq >= 4.
    Calls info_polymorphism().

        Parameters:
            alignment (list): list of bases
        
        Output :
            list: positions of informative polymorphisms
    """
    polym_positions = []

    for i in range(len(alignment)):
        if info_polymorphism(alignment[i]):
            polym_positions.append(i)

    return polym_positions


def four_gametes_test(bases1 :list, bases2 :list):
    """
    Whether 2 loci are incompatible by computing a 4-gametes
    test on 2 informative polymorphisms.

        Parameters:
            bases1 (list)
            bases2 (list)
        
        Output:
            boolean: the 2 loci are incompatible
    """
    assert len(bases1) == len(bases2), "Not the same number of elements."

    c = []
    for i in range(len(bases1)):
        c.append((bases1[i], bases2[i]))

    return(len(set(c)) > 3)


def incompatible_polym(alignment :list,
                       polym_pos :list):
    """
    From an alignment, returns a list of pairwise incompatible polymorphisms.
    Calls polym_positions() and four_gametes_test().

        Parameters:
            alignment (list of lists)
            polym_pos (list): polym_positions() output  
        Output:
            list of tuples: positions of incompatiblities    
    """
    incomp_polym = []
    #polym_pos = polym_positions(alignment)

    for i in range(len(polym_pos)):
        first_ind = polym_pos[i]
        for j in range(i, len(polym_pos)):
            second_ind = polym_pos[j]
            if four_gametes_test(alignment[first_ind], alignment[second_ind]):
                incomp_polym.append((first_ind, second_ind))

    return incomp_polym


def chop_list(list_pos :list):
    """
    Seeks the shortest interval that is sufficient to explain the
    incompatibilities.

    Paremeters:
        list_pos (list of tuples): pairwise positions of incompatible loci
    
    Output:
        list : chopped list

    ! : "We retrieve all intervals and sort them in increasing order of site
    positions along the genome" -> DONE by constructing, not verified anywhere
    """
    curr = list_pos[-1]

    for interval in list_pos[::-1][1:]:
        # the current interval is contained in the other one
        if curr[0] >= interval[0] and curr[0] < interval[1] and \
            curr[1] <= interval[1] :
            list_pos.remove(interval)    

        # overlap
        elif curr[0] >= interval[0] and curr[0] < interval[1] and \
            curr[1] >= interval[1]:
            list_pos.remove(curr)
            list_pos.remove(interval)
            new_interv = (curr[0], interval[1])
            curr = new_interv
            list_pos.append(new_interv)

        else : # move on
            curr = interval

    return(list_pos)


def place_mld_breakpoints(incompt_pos :list):
    """
    Adds an MLD breakpoint in the middle of an incompatible interval.

        Parameters:
            incompt_pos (list): chop_list() output
    
        Output:
            list: indexes of MLD breakpoints in the sequence
    """
    return [int(inter[0] + (inter[1] - inter[0])/2) for inter in incompt_pos]


def mld_length_distrib(mld_pos :list):
    """
    Plots an histogram of MLD bloc lengths.

        Parameters:
            mld_pos (list): mld positions,
                typically place_mld_breakpoints() output
    
    Requires matplotlib.pyplot.
    """
    #TODO
