from generatetags.find_possible import get_tag_interval, tagint

def test_tag_interval_primer_is_left():
    p = tagint('chr1', 0, 10, None, '+')
    re = tagint('chr1', 15, 20, None, '+')

    expected = tagint('chr1', 0, 20, None, '+')
    res = get_tag_interval(p, re)
    assert expected == res

def test_tag_interval_primer_is_left_and_overlaps():
    p = tagint('chr1', 0, 10, None, '+')
    re = tagint('chr1', 9, 20, None, '+')

    expected = tagint('chr1', 0, 20, None, '+')
    res = get_tag_interval(p, re)
    assert expected == res


def test_tag_interval_primer_is_right():
    p = tagint('chr21', 41458235, 41458255, 'BACE2_rev', '-')
    re = tagint('chr21', 41458227, 41458231, 'rsite', '+')

    expected = tagint('chr21', 41458227, 41458255, 'BACE2_rev', '-')
    res = get_tag_interval(p, re, name='BACE2_rev')
    print expected
    print res
    assert expected == res

def test_tag_interval_primer_is_right_and_overlaps():
    p = tagint('chr1', 9, 20, None, '-')
    re = tagint('chr1', 0, 10, None, '+')

    expected = tagint('chr1', 0, 20, None, '-')
    res = get_tag_interval(p, re)
    print expected
    print res
    assert expected == res
