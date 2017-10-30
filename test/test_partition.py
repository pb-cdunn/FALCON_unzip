import falcon_unzip.mains.bam_partition_and_merge as M
import StringIO
import collections
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
