import os
import falcon_unzip.io as M

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
