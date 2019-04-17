import falcon_unzip.mains.db_to_ccs_id as mod
import pytest
import helpers
import os
import filecmp


def test_help():
    try:
        mod.main(['prog', '--help'])
    except SystemExit:
        pass

def test_main_1(request):

    test_data = os.path.join(helpers.get_test_data_dir(), '0-phasing/')

    argv = ['prog',
            '--lookup'      , test_data + 'readname_lookup.txt',
            '--rid-to-phase', test_data + 'rid_to_phase.tmp',
            '--rid-to-ctg'  , test_data + 'rid_to_ctg.txt',
            '--ctg'         , '000000F',
            '--output'      , 'test_db_to_ccs_id_res.txt'
            ]

    mod.main(argv)

    assert(filecmp.cmp(test_data + 'rid_to_phase', 'test_db_to_ccs_id_res.txt'))
