import pytest
import falcon_unzip.mains.bam_partition_and_merge as M
import StringIO
import collections
import logging
import os

read2ctg = {
    'r11': 'c1',
    'r12': 'c1',
    'r21': 'c2',
    'r31': 'c3',
    'r32': 'c3',
    'r33': 'c3',
}


def test_partition_ctgs_1():
    groups = M.partition_ctgs(read2ctg, 1)
    assert groups == [set(['c1', 'c2', 'c3'])]


def test_partition_ctgs_2():
    groups = M.partition_ctgs(read2ctg, 2)
    assert groups == [set(['c1', 'c2']), set(['c3'])]


def test_partition_ctgs_3():
    groups = M.partition_ctgs(read2ctg, 3)
    assert groups == [set(['c3']), set(['c1']), set(['c2'])]
    # Order matters. They tend to be sorted by number of reads.


def test_partition_ctgs_4():
    groups3 = M.partition_ctgs(read2ctg, 3)
    groups4 = M.partition_ctgs(read2ctg, 4)
    assert groups3 == groups4


def test_get_zmw():
    assert M.get_zmw('foo/123/anything') == 'foo/123'
    assert M.get_zmw('foo/123') == 'foo/123'
    with pytest.raises(Exception) as e:
        M.get_zmw('foo')
    with pytest.raises(Exception) as e:
        M.get_zmw('foo/123/anything/extra')

def test_get_zmw2ctg(caplog):
    caplog.set_level(logging.WARN)
    read2ctg = {
        'foo/123/0_1': 'ABC',
        'foo/123/2_3': 'ABC',
        'foo/123/4_5': 'DEF',
        'bar/321/0_1': 'ABC',
        'bar/654/0_1': 'DEF',
    }
    zmw2ctg = M.get_zmw2ctg(read2ctg)
    assert zmw2ctg['foo/123'] in ('ABC', 'DEF')
    assert zmw2ctg['bar/321'] == 'ABC'
    assert zmw2ctg['bar/654'] == 'DEF'
    assert 'foo/123' in caplog.text
    assert 'bar' not in caplog.text
