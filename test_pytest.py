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
    # TO DO
    pass                


def test_find_the_culprit():
    bases1 = ['A','A','a', 'a', 'A', 'A', 'a', 'A', 'a']
    bases2 = ['B','b','b', 'b', 'B', 'b', 'b', 'B', 'B']
    # ('A', 'B') : 3, ('A','b') : 2, ('a','b') : 2, ('B', 'a') : 1
    # -> 5
    assert find_the_culprit(bases1, bases2) == 8
