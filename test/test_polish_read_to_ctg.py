import falcon_unzip.mains.polish_read_to_ctg as mod
import pytest
import os
import filecmp


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
    out_fn = str(tmpdir.join('read2ctg.txt'))
    def td1(fn):
        return os.path.join(request.fspath.dirname, '..', 'test_data', '0-phasing', fn)
    def td2(fn):
        return os.path.join(request.fspath.dirname, '..', 'test_data', '4-polishing', fn)
    argv = ['prog',
            '--lookup-fn'      , td1('readname_lookup.txt'),
            '--rid-to-phase-fn', td2('rid_to_phase.all'),
            '--edges-fn'       , td2('combined_edges.txt'),
            '--out-read-to-ctg', out_fn
            ]

    mod.main(argv)

    filecheck(td2('read2ctg.txt'), out_fn)
