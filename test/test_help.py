import importlib
import pytest


@pytest.mark.parametrize('mod', [
    'bam_partition_and_merge',
    'bam_segregate',
    'dedup_h_tigs',
    'get_read2ctg',
    'get_read_hctg_map',
    'graphs_to_h_tigs',
    'ovlp_filter_with_phase',
    'phased_ovlp_to_graph',
    'phasing',
    'phasing_readmap',
    'rr_hctg_track',
    'start_quiver',
    'start_unzip',
    'unzip_gen_gfa_v1',

    # These are not actually used, but we can still check them.
    'get_ctg2bam_map',
    'select_reads_from_bam',
    'track_reads_htigs0',
])
def test(mod):
    module = importlib.import_module('falcon_unzip.mains.{}'.format(mod))
    with pytest.raises(SystemExit) as excinfo:
        module.main(['prog', '--help'])
    assert 0 == excinfo.value.code
