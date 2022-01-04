from four_gametes import *

def test_info_polymorphism():
    assert info_polymorphism(['A', 'A', 'a', 'A', 'a']) == True
    assert info_polymorphism(['B', 'B', 'B', 'b', 'B']) == False


def test_polym_positions():
    assert polym_positions(['ABCDEF',
                            'aBcDeF',
                            'abcdeF',
                            'ABcdEF',
                            'ABCdEF']) == [0, 2, 3, 4]


def test_four_gametes_test():
    assert four_gametes_test(['A','A','a', 'a', 'A', 'A'],
                             ['B','b','b', 'b', 'B', 'b']) == False
    assert four_gametes_test(['A','A','a', 'a', 'A', 'A','a'],
                             ['B','b','b', 'b', 'B', 'b', 'B']) == True


def test_incompatible_polym():
    # TO DO
    pass


def test_chop_list():
    assert chop_list([[1,4], [2,3]]) == [[2,3]]
    assert chop_list([(1,3), (5, 20), (6,12)]) == [(1,3), (6,12)]
    assert chop_list([(1,3), (5,17), (5,20), (6, 12)]) == [(1,3), (6,12)]
    assert chop_list([(3,10), (6, 12)]) == [(6,10)]
    assert chop_list([(3, 8), (3, 10), (5, 17), (5, 20), (6, 12)]) == [(6, 8)]
    assert chop_list([(2,3), (2, 6), (3, 5), (3, 10), (5, 17), (5, 20), \
        (6, 12)]) == [(2, 3), (3, 5), (6, 10)]
    assert chop_list([(1,2), (1,4), (2,3), (2, 6), (3, 5), (3, 10), (5, 17), \
         (5, 20), (6, 12)]) == [(1,2), (2, 3), (3, 5), (6, 10)]



def test_find_the_culprit():
    bases1 = ['A','A','a', 'a', 'A', 'A', 'a', 'A', 'a']
    bases2 = ['B','b','b', 'b', 'B', 'b', 'b', 'B', 'B']
    # ('A', 'B') : 3, ('A','b') : 2, ('a','b') : 2, ('B', 'a') : 1
    # -> 8
    assert find_the_culprit(bases1, bases2) == 8




