import os
import pytest
import subprocess
import falcon_unzip.io as M

import logging
#logging.basicConfig(level=logging.INFO)

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


def test_update_env_from_config(tmpdir, monkeypatch):
    which_cmd = M.capture('which which').strip()

    monkeypatch.setenv('PATH', '')
    config = {}
    M.update_env_from_config(config, None)

    monkeypatch.setenv('PATH', os.path.dirname(which_cmd))
    with pytest.raises(Exception) as excinfo:
        config = {}
        M.update_env_from_config(config, None)
    assert "Call 'which " in str(excinfo.value)

    cmds = [
        'blasr', 'samtools', 'pbalign', 'variantCaller',
        'nucmer', 'show-coords',
        'fc_rr_hctg_track2.exe',
    ]
    for cmd in cmds:
        p = tmpdir.join(cmd)
        p.write('#!/bin/bash')
        p.chmod(0o777)

    config = {
        'smrt_bin': str(tmpdir)
    }
    M.update_env_from_config(config, None)
