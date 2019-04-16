import falcon_unzip.mains.db_to_ccs_id as mod
import pytest
import filecmp

def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass

def test_main_1(request):
    print("dir:", request.fspath.dirname)
    argv = ['prog',
            '--lookup'      , '../test_data/0-phasing/readname_lookup.txt',
            '--rid-to-phase', '../test_data/0-phasing/rid_to_phase.tmp',
            '--rid-to-ctg'  , '../test_data/0-phasing/rid_to_ctg.txt',
            '--ctg'         , '000000F',
            '--output'      , 'test_db_to_ccs_id_res.txt'
            ]

    mod.main(argv)

    assert(filecmp.cmp('../test_data/0-phasing/rid_to_phase', 'test_db_to_ccs_id_res.txt'))
