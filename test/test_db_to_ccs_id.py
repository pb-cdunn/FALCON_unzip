import falcon_unzip.mains.db_to_ccs_id as mod
import pytest
import os


def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass

def filecheck(got_fn, expected_fn):
    got = open(got_fn).read()
    expected = open(expected_fn).read()
    assert got == expected

def test_main_1(request, tmpdir):
    out_fn = str(tmpdir.join('test_db_to_ccs_id_res.txt'))
    def td(fn):
        return os.path.join(request.fspath.dirname, '..', 'test_data', '0-phasing', fn)
    argv = ['prog',
            '--lookup'      , td('readname_lookup.txt'),
            '--rid-to-phase', td('rid_to_phase.tmp'),
            '--rid-to-ctg'  , td('rid_to_ctg.txt'),
            '--ctg'         , '000000F',
            '--output'      , out_fn,
            ]

    mod.main(argv)

    filecheck(td('rid_to_phase'), out_fn)
