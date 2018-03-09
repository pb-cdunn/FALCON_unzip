import commands
import json
import logging
import msgpack
import os
import pprint
import re
import subprocess

try:
    # pylint: disable=no-name-in-module, import-error, fixme, line-too-long
    from pysam.calignmentfile import AlignmentFile
except ImportError:
    # pylint: disable=no-name-in-module, import-error, fixme, line-too-long
    from pysam.libcalignmentfile import AlignmentFile

LOG = logging.getLogger()


def log(*msgs):
    LOG.info(' '.join(repr(m) for m in msgs))


def validate_samtools(samtools_output):
    """Given the result of a bare call to 'samtools',
    prove our $PATH is using at least version 1.3.0
    """
    re_Version = re.compile(r'Version: (?P<major>\d+)\.(?P<minor>\d+)\.(?P<subminor>\d+)')
    mo = re_Version.search(samtools_output)
    if not mo:
        msg = samtools_output + '\n---\nCould not discern version of samtools. We need at least 1.3.0. Good luck!'
        LOG.warning(msg)
        return
    major = int(mo.group('major'))
    minor = int(mo.group('minor'))
    if major < 1 or (major == 1 and minor < 3):
        msg = 'samtools is version {}.{}.{}'.format(
                major, minor, mo.group('subminor'))
        msg += ', but we require >= 1.3.0'
        raise Exception(msg)


def validate_config(config):
    # This simple and quick check catches common problems early.
    # This code might go somewhere else someday.
    smrt_bin_cmds = [
        'blasr', 'samtools', 'pbalign', 'variantCaller',
    ]
    path_cmds = [
        'nucmer',
        'show-coords',
        'fc_rr_hctg_track2.exe',
    ]
    LOG.info('PATH={}'.format(os.environ['PATH']))
    try:
        capture('which which')
    except Exception:
        LOG.warning('Could not find "which" command. Skipping checks for "blasr", etc.')
        return
    for cmd in smrt_bin_cmds + path_cmds:
        syscall('which ' + cmd)
    syscall('show-coords -h')
    syscall('nucmer --version')

    samtools_output = commands.getoutput('samtools')
    validate_samtools(samtools_output)


def update_env_from_config(config, fn):
    LOG.info('From {!r}, config={}'.format(fn, pprint.pformat(config)))
    smrt_bin = config.get('smrt_bin')
    if smrt_bin:
        PATH = '{}:{}'.format(os.environ['PATH'], smrt_bin)
        os.environ['PATH'] = PATH
    validate_config(config)


def mkdirs(*dirnames):
    for dirname in dirnames:
        if not dirname:
            continue # '' => curdir
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
            if len(dirnames) == 1:
                log('mkdir -p {!r}'.format(dirnames[0]))


def eng(number):
    return '{:.1f}MB'.format(number / float(2**20))


def read_as_msgpack(stream):
    content = stream.read()
    log('  Read {} as msgpack'.format(eng(len(content))))
    return msgpack.loads(content)


def read_as_json(stream):
    content = stream.read()
    log('  Read {} as json'.format(eng(len(content))))
    return json.loads(content)


def write_as_msgpack(stream, val):
    content = msgpack.dumps(val)
    log('  Serialized to {} as msgpack'.format(eng(len(content))))
    stream.write(content)


def write_as_json(stream, val):
    content = json.dumps(val, indent=2, separators=(',', ': '))
    log('  Serialized to {} as json'.format(eng(len(content))))
    stream.write(content)


def deserialize(fn):
    log('Deserializing from {!r}'.format(fn))
    with open(fn) as ifs:
        log('  Opened for read: {!r}'.format(fn))
        if fn.endswith('.msgpack'):
            val = read_as_msgpack(ifs)
        elif fn.endswith('.json'):
            val = read_as_json(ifs)
        else:
            raise Exception('Unknown extension for {!r}'.format(fn))
    log('  Deserialized {} records'.format(len(val)))
    return val


def serialize(fn, val):
    """Assume dirname exists.
    """
    log('Serializing {} records'.format(len(val)))
    mkdirs(os.path.dirname(fn))
    with open(fn, 'w') as ofs:
        log('  Opened for write: {!r}'.format(fn))
        if fn.endswith('.msgpack'):
            write_as_msgpack(ofs, val)
        elif fn.endswith('.json'):
            write_as_json(ofs, val)
            ofs.write('\n') # for vim
        else:
            raise Exception('Unknown extension for {!r}'.format(fn))


def yield_bam_fn(input_bam_fofn_fn):
    log('Reading BAM names from FOFN {!r}'.format(input_bam_fofn_fn))
    fofn_basedir = os.path.normpath(os.path.dirname(input_bam_fofn_fn))

    def abs_fn(maybe_rel_fn):
        if os.path.isabs(maybe_rel_fn):
            return maybe_rel_fn
        else:
            return os.path.join(fofn_basedir, maybe_rel_fn)
    for row in open(input_bam_fofn_fn):
        yield abs_fn(row.strip())


def yield_abspath_from_fofn(fofn_fn):
    """Yield each filename.
    Relative paths are resolved from the FOFN directory.
    """
    try:
        fns = deserialize(fofn_fn)
    except:
        #LOG('las fofn {!r} does not seem to be JSON; try to switch, so we can detect truncated files.'.format(fofn_fn))
        fns = open(fofn_fn).read().strip().split()
    try:
        basedir = os.path.dirname(fofn_fn)
        for fn in fns:
            if not os.path.isabs(fn):
                fn = os.path.abspath(os.path.join(basedir, fn))
            yield fn
    except Exception:
        LOG.error('Problem resolving paths in FOFN {!r}'.format(fofn_fn))
        raise


def syscall(call, nocheck=False):
    """Raise Exception in error, unless nocheck==True
    """
    LOG.info('$(%s)' % repr(call))
    rc = os.system(call)
    msg = 'Call %r returned %d.' % (call, rc)
    if rc:
        LOG.warning(msg)
        if not nocheck:
            raise Exception(msg)
    else:
        LOG.debug(msg)
    return rc


def capture(cmd):
    """Return stdout, fully captured.
    Wait for subproc to finish.
    Raise if empty.
    Raise on non-zero exit-code.
    """
    LOG.info('$ {} >'.format(cmd))
    output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    if not output:
        msg = '{!r} failed to produce any output.'.format(cmd)
        LOG.warning(msg)
    return output


def rm(*f):
    syscall('rm -f {}'.format(' '.join(f)))


def touch(f):
    syscall('touch {}'.format(f))


def filesize(fn):
    return os.stat(fn).st_size


def exists_and_not_empty(fn):
    if not os.path.exists(fn):
        return False
    if 0 == filesize(fn):
        LOG.debug('File {} is empty.'.format(fn))
        return False
    return True


def substitute(yourdict):
    """
    >>> list(sorted(substitute({'a': '_{b}_', 'b': 'X'}).items()))
    [('a', '_X_'), ('b', 'X')]
    """
    mydict = dict(yourdict)
    for k, v in yourdict.items():
        if '{' in v:
            mydict[k] = v.format(**mydict)
    return mydict
