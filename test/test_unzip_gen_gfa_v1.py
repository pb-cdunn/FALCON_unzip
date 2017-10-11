import falcon_unzip.unzip_gen_gfa_v1 as mod
import helpers
import pytest
import os


def test_main_1(capsys):
    """
    Tests writing a GFA from both sg.gexf and tiling paths.
    Outputs read and contig sequences.
    """
    test_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')

    argv = ['prog',
            '--preads-fasta', os.path.join(test_dir, 'preads4falcon.fasta'),
            '--p-ctg-fasta', os.path.join(test_dir, 'p_ctg.fa'),
            '--h-ctg-fasta', os.path.join(test_dir, 'h_ctg.fa'),
            '--unzip-root', os.path.join(test_dir, './'),
            # '--tiling',
            '--write-reads',
            '--write-contigs',
            '--min-p-len', '0',
            '--min-h-len', '0',
            ]
    mod.main(argv)
    out, err = capsys.readouterr()
    expected_path = os.path.join(test_dir, 'expected-1-sg-r-c.gfa')
    expected = open(expected_path).read()
    assert(out == expected)


def test_main_2(capsys):
    """
    Tests writing a GFA only from tiling paths.
    Outputs read and contig sequences.
    """
    test_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')

    argv = ['prog',
            '--preads-fasta', os.path.join(test_dir, 'preads4falcon.fasta'),
            '--p-ctg-fasta', os.path.join(test_dir, 'p_ctg.fa'),
            '--h-ctg-fasta', os.path.join(test_dir, 'h_ctg.fa'),
            '--unzip-root', os.path.join(test_dir, './'),
            '--tiling',
            '--write-reads',
            '--write-contigs',
            '--min-p-len', '0',
            '--min-h-len', '0',
            ]
    mod.main(argv)
    out, err = capsys.readouterr()
    expected_path = os.path.join(test_dir, 'expected-2-tiling-r-c.gfa')
    expected = open(expected_path).read()
    assert(out == expected)


def test_main_3(capsys):
    """
    Tests writing a GFA only from tiling paths.
    Outputs contig sequences but not reads.
    """
    test_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')

    argv = ['prog',
            '--preads-fasta', os.path.join(test_dir, 'preads4falcon.fasta'),
            '--p-ctg-fasta', os.path.join(test_dir, 'p_ctg.fa'),
            '--h-ctg-fasta', os.path.join(test_dir, 'h_ctg.fa'),
            '--unzip-root', os.path.join(test_dir, './'),
            '--tiling',
            # '--write-reads',
            '--write-contigs',
            '--min-p-len', '0',
            '--min-h-len', '0',
            ]
    mod.main(argv)
    out, err = capsys.readouterr()
    expected_path = os.path.join(test_dir, 'expected-3-tiling-no_r-c.gfa')
    expected = open(expected_path).read()
    assert(out == expected)


def test_main_4(capsys):
    """
    Tests writing a GFA only from tiling paths.
    Does not output contig or read sequences.
    """
    test_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')

    argv = ['prog',
            '--preads-fasta', os.path.join(test_dir, 'preads4falcon.fasta'),
            '--p-ctg-fasta', os.path.join(test_dir, 'p_ctg.fa'),
            '--h-ctg-fasta', os.path.join(test_dir, 'h_ctg.fa'),
            '--unzip-root', os.path.join(test_dir, './'),
            '--tiling',
            # '--write-reads',
            # '--write-contigs',
            '--min-p-len', '0',
            '--min-h-len', '0',
            ]
    mod.main(argv)
    out, err = capsys.readouterr()
    expected_path = os.path.join(test_dir, 'expected-4-tiling-no_r-no_c.gfa')
    expected = open(expected_path).read()
    assert(out == expected)


def test_main_5(capsys):
    """
    Tests writing a GFA from both sg.gexf and tiling paths.
    Does not output contig or read sequences.
    """
    test_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')

    argv = ['prog',
            '--preads-fasta', os.path.join(test_dir, 'preads4falcon.fasta'),
            '--p-ctg-fasta', os.path.join(test_dir, 'p_ctg.fa'),
            '--h-ctg-fasta', os.path.join(test_dir, 'h_ctg.fa'),
            '--unzip-root', os.path.join(test_dir, './'),
            # '--tiling',
            # '--write-reads',
            # '--write-contigs',
            '--min-p-len', '0',
            '--min-h-len', '0',
            ]
    mod.main(argv)
    out, err = capsys.readouterr()
    expected_path = os.path.join(test_dir, 'expected-5-sg-no_r-no_c.gfa')
    expected = open(expected_path).read()
    assert(out == expected)


def test_main_6(capsys):
    """
    Tests writing a GFA only from tiling paths.
    Does not output contig or read sequences.
    Filters tiling paths according to specified lengths.
    """
    test_dir = os.path.join(helpers.get_test_data_dir(), 'gfa-1')

    argv = ['prog',
            '--preads-fasta', os.path.join(test_dir, 'preads4falcon.fasta'),
            '--p-ctg-fasta', os.path.join(test_dir, 'p_ctg.fa'),
            '--h-ctg-fasta', os.path.join(test_dir, 'h_ctg.fa'),
            '--unzip-root', os.path.join(test_dir, './'),
            '--tiling',
            # '--write-reads',
            # '--write-contigs',
            '--min-p-len', '10000',
            '--min-h-len', '10000',
            ]
    mod.main(argv)
    out, err = capsys.readouterr()
    expected_path = os.path.join(test_dir, 'expected-6-tiling-no_r-no_c-minlen.gfa')
    expected = open(expected_path).read()
    assert(out == expected)
