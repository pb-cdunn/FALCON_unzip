import falcon_unzip.mains.polish_unphased_readmapping as mod
import pytest
import os

def filecheck(got_fn, expected_fn):
    got = open(got_fn).read()
    expected = open(expected_fn).read()
    assert got == expected

def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass

def test_main_1(request, tmpdir):
    out_rn = str(tmpdir.join('read_names.txt'))
    out_sn = str(tmpdir.join('ref_names.txt'))
    def td(fn):
        return os.path.join(request.fspath.dirname, '..', 'test_data', '4-polishing', fn)
    argv = ['prog',
            '--ctg'           , '000000F',
            '--fai'           , td('combined_ph.fa.fai'),
            '--out-read-names', out_rn, 
            '--out-ref-names' , out_sn,
            '--read-to-ctg'   , td('read2ctg2.txt')
            ]

    mod.main(argv)

    filecheck(td('read_names.txt'), out_rn)
    filecheck(td('ref_names.txt'), out_sn)
