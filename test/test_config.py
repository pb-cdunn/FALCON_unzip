from falcon_unzip import unzip as M
from falcon_unzip import io
import pytest

def test_validate_input_bam_fofn(tmpdir):
    with io.cd(str(tmpdir)):
        cfg = {
            'Unzip': {
            }
        }
        with pytest.raises(AssertionError) as exc:
            M.validate_input_bam_fofn(cfg, 'dummy')
        assert 'You must provide "input_bam_fofn"' in str(exc.value)


        cfg = {
            'Unzip': {
                'input_bam_fofn': 'foo.fofn'
            }
        }
        with pytest.raises(Exception) as exc:
            M.validate_input_bam_fofn(cfg, 'dummy')
        assert "No such file or directory: 'foo.fofn'" in str(exc.value)

        io.touch('foo.fofn')
        M.validate_input_bam_fofn(cfg, 'dummy')

        with open('foo.fofn', 'w') as sout:
            sout.write("""\
bar.txt
baz.txt
""")
        with pytest.raises(Exception) as exc:
            M.validate_input_bam_fofn(cfg, 'dummy')
        assert 'All filenames in input_bam FOFN ("foo.fofn") must end in ".bam".' in str(exc.value)

        with open('foo.fofn', 'w') as sout:
            sout.write("""\
bar.bam
baz.bam
""")
        M.validate_input_bam_fofn(cfg, 'dummy')
        # It does not actually validate existence.
