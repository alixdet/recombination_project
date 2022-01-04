"""
Long project
SMILE unit

Ref : 
- Kerdoncuff, Lambert & Achaz, 2020
"""

from collections import Counter

def info_polymorphism(bases :list):
    """
    From a list of at least 4 elements, returns whether the
    polymorphism is informative
    i.e. 2 alleles at least twice

    ex, on a given locus :
        - A A a A a -> True
        - B B B B b -> False
    
    ! : consider case of random mutation on informative polymorphism
        when # seq >= 5

    if remains that simple, delete function and insert code instead
    """
    assert len(bases) >= 4, "4 or more elments required"

    # maybe find a way w/o using Counter
    values = Counter(bases).values()

    return len(values) >= 2 and 1 not in values


def polym_positions(alignment :list):
    len_alignment = len(alignment[0])
    for seq in alignment[1:]:
        assert len(seq) == len_alignment, "Sequences of different sizes"
    
    polym_positions = []

    for i in range(len_alignment):
        if info_polymorphism([seq[i] for seq in alignment]):
            polym_positions.append(i)

    return polym_positions


def four_gametes_test(bases1 :list, bases2 :list):
    """
    4-gametes test on 2 informative polymorphisms loci

    Returns True if 2 loci are incompatible, else False
    """
    assert len(bases1) == len(bases2)

    c = []
    for i in range(len(bases1)):
        c.append((bases1[i], bases2[i]))

    return(len(set(c)) > 3)


def incompatible_polym(alignment :list):
    """
    Return list of pairwise incompatible polymorphisms
    """
    incomp_polym = []
    polym_pos = polym_positions(alignment)

    for i in range(len(polym_pos)):
        for j in range(len(polym_pos[i:])):
            if four_gametes_test([seq[i] for seq in alignment],
                                [seq[j] for seq in alignment]):
                incomp_polym.append((i,j))

    return incomp_polym

    """tests
    TO DO
        """

def chop_list(list_pos :list):
    """
    ...
    """


def find_the_culprit(bases1 :list, bases2 :list):
    """
    based on four_gametes_test()

    with 2 incomaptible polymorphism, removes the positions (i.e. the sequence)
    that are responsible for the incompatibility (i.e. whom combination appears
    the least)

    returns indexes of compatible sequences
    """
    #assert len(bases1) == len(bases2)
    comb = {}
    for i in range(len(bases1)):
        comb.setdefault((bases1[i], bases2[i]), [])
        comb[(bases1[i], bases2[i])].append(i)

    if len(set(comb)) == 4 :
        min_occ = min([len(values) for values in comb.values()])
        min_values = [v for v in comb.values() if len(v) == min_occ][0]
        return min_values[0]
        
    return False # not sure what to return


def longuest_compat_chunk_notime(tot_alignment :list, start_pos :int):
    """
    DOES NOT CONSIDER THE TIME

    Resize each time sequences are incompatible

    Counts the # of mutations

    Starts at a polymorphic locus
    """
    bases0 = [seq[start_pos] for seq in tot_alignment]
    mutations = 0

    len_alignment = len(tot_alignment[0])
    for seq in tot_alignment[1:]:
        assert len(seq) == len_alignment, "Sequences of different sizes"
    
    alignment = tot_alignment # we start with the entire alignment

    i = start_pos + 1
    while (i < len_alignment and len(alignment) > 3):
        # ! : len(alignment) = # of seq =/= length of sequences aligned
        bases = [seq[i] for seq in alignment]
        # is the locus polymorphic ?
        # !! : do 2 mutations on same locus count as 2 mutation events ?
        # we start by saying no
        if len(Counter(bases) > 1):
            mutations += 1

            # terrible way, culprit is or a list or False !
            culprit = find_the_culprit(bases0, bases)
            if culprit:
                del alignment[culprit]

        i += 1
    
    return (start_pos, i)
    ### TO TEST !


def longuest_compat_recent_bloc(tot_alignment :list, start_pos :int):
    """
    Resize each time sequences are incompatible

    Counts the # of mutations

    Starts at a polymorphic locus
    """
    bases0 = [seq[start_pos] for seq in tot_alignment]
    mutations = 0
    end_time = 999 #find a format, value, ...

    len_alignment = len(tot_alignment[0])
    for seq in tot_alignment[1:]:
        assert len(seq) == len_alignment, "Sequences of different sizes"
    
    alignment = tot_alignment # we start with the entire alignment

    i = start_pos + 1
    time = 0

    while (i < len_alignment and len(alignment) > 3 and time < end_time):
        # ! : len(alignment) = # of seq =/= length of sequences aligned
        bases = [seq[i] for seq in alignment]
        # is the locus polymorphic ?
        # !! : do 2 mutations on same locus count as 2 mutation events ?
        # we start by saying no
        if len(Counter(bases) > 1):
            mutations += 1

            # terrible way, culprit is or a list or False !
            culprit = find_the_culprit(bases0, bases)
            if culprit:
                del alignment[culprit]

        #time = mutations x mutation_rate x #seq
        i += 1
    
    return (start_pos, i)