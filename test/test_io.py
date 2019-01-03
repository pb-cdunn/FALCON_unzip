from StringIO import StringIO
import os
import pytest
import subprocess
import falcon_unzip.io as M

import logging
# logging.basicConfig(level=logging.INFO)


def test_exists_and_not_empty(tmpdir):
    with tmpdir.as_cwd():
        foo = 'foo'
        assert not os.path.exists(foo)
        assert not M.exists_and_not_empty(foo)
        M.touch(foo)
        assert os.path.exists(foo)
        assert not M.exists_and_not_empty(foo)
        with open(foo, 'w') as stream:
            stream.write('ABC')
        assert 3 == M.filesize(foo)
        assert M.exists_and_not_empty(foo)
        M.rm(foo)
        assert not M.exists_and_not_empty(foo)


def test_eng():
    assert '0.1MB' == M.eng(100000)
    assert '0.1MB' == M.eng(2**17)


def test_validate_config(tmpdir, monkeypatch):
    #which_cmd = M.capture('which which').strip()
    cmds = [
        'blasr', 'samtools', 'pbalign', 'variantCaller',
        'minimap2',
        'nucmer', 'show-coords',
        'fc_rr_hctg_track2.exe',
    ]
    for cmd in cmds:
        p = tmpdir.join(cmd)
        p.write('#!/bin/bash')
        p.chmod(0o777)
    config = dict()
    monkeypatch.setenv('PATH', str(tmpdir))
    M.validate_config(config)


def test_validate_samtools():
    with pytest.raises(Exception) as excinfo:
        M.validate_samtools('Version: 1.2.11')
    assert 'but we require' in str(excinfo.value)

    M.validate_samtools('Version: 1.3.1')
    M.validate_samtools('Version: 2.1.0')
    M.validate_samtools('Version: 1.6 (using htslib 1.7)')


# In the future, this can be used to learn the directory of the current file.
@pytest.fixture
def thisdir(request):
    return os.path.dirname(request.module.__file__)

fofn = """\
/pbi/1.subreads.bam
/pbi/2.subreads.bam

"""

def test_yield_bam_fn(tmpdir):
    assert fofn.endswith('\n\n')
    with tmpdir.as_cwd():
        fofn_fn = './foo.fn'
        with open(fofn_fn, 'w') as stream:
            stream.write(fofn)
        fns = list(M.yield_bam_fn(fofn_fn))
        assert fns == ['/pbi/1.subreads.bam', '/pbi/2.subreads.bam']

    with tmpdir.as_cwd():
        bad_fofn = fofn + "/foo/bar.bai"
        fofn_fn = './foo.fn'
        with open(fofn_fn, 'w') as stream:
            stream.write(bad_fofn)
        with pytest.raises(Exception) as excinfo:
            list(M.yield_bam_fn(fofn_fn))
        assert 'All filenames in input_bam FOFN' in str(excinfo.value)
        assert '/foo/bar.bai' in str(excinfo.value)
        assert '/foo.fn' in str(excinfo.value)
