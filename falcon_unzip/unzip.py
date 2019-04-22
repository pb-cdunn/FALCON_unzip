from falcon_kit import run_support as support
from falcon_kit import snakemake
from falcon_kit.FastaReader import FastaReader
from pypeflow.simple_pwatcher_bridge import (
    PypeProcWatcherWorkflow,
)
from .tasks import unzip as tasks_unzip, ccs as tasks_ccs
from . import io
import glob
import logging
import os
import re
import time
import falcon_kit.run_support
import falcon_kit.functional

LOG = logging.getLogger(__name__)


def unzip_all(config, unzip_config_fn):
    job_defaults = config['job.defaults']
    use_tmpdir = job_defaults['use_tmpdir'] # None/False is fine.
    wf = PypeProcWatcherWorkflow(
        job_defaults=job_defaults,
        use_tmpdir=use_tmpdir,
    )
    tasks_unzip.run_workflow(wf, config, unzip_config_fn)

def ccs_all(config, unzip_config_fn):
    job_defaults = config['job.defaults']
    use_tmpdir = job_defaults['use_tmpdir'] # None/False is fine.
    wf = PypeProcWatcherWorkflow(
        job_defaults=job_defaults,
        use_tmpdir=use_tmpdir,
    )
    tasks_ccs.run_workflow(wf, config, unzip_config_fn)

def update_config_from_sections(config, cfg):
    allowed_sections = set([
            'General',
            'Unzip',
            'job.step.unzip.track_reads',
            'job.step.unzip.phasing',
            'job.step.unzip.blasr_aln',
            'job.step.unzip.hasm',
            'job.step.unzip.quiver',
            'job.defaults',
    ])
    all_sections = set(k for k,v in cfg.items() if isinstance(v, dict))
    unexpected = all_sections - allowed_sections
    if unexpected:
        msg = 'You have {} unexpected cfg sections: {}'.format(
            len(unexpected), unexpected)
        raise Exception(msg)
    config.update(cfg)
    # Guarantee they all exist.
    for sec in allowed_sections:
        if sec not in config:
            config[sec] = dict()

def update_defaults(config):
    """Update dict with any missing defaults.
    """
    def set_default(section, key, val):
        cfg = config[section]
        if key not in cfg:
            cfg[key] = val
    set_default('Unzip', 'polish_vc_ignore_error', False)
    set_default('Unzip', 'polish_use_blasr', False)
    set_default('Unzip', 'polish_include_zmw_all_subreads', False)

    # Fix up known boolean config-values, which could be strings.
    for section, bool_key in (
            ('Unzip', 'polish_vc_ignore_error'),
            ('Unzip', 'polish_use_blasr'),
            ('Unzip', 'polish_include_zmw_all_subreads'),
            ):
        cfg = config[section]
        cfg[bool_key] = falcon_kit.functional.cfg_tobool(cfg[bool_key])

def parse_cfg_file(config_fn):
    """Return as dict.
    This overlooks some ad-hoc config.
    """
    # New: Parse sections (case-sensitively), into sub-dicts.
    config = dict()
    with open(config_fn) as stream:
        cfg2 = falcon_kit.run_support.parse_cfg_with_sections(stream)
        update_config_from_sections(config, cfg2)
    falcon_kit.run_support.update_job_defaults_section(config)
    update_defaults(config)
    return config

def parse_config(config_fn):
    config = parse_cfg_file(config_fn)

    # We might need these in the top-level.
    config.setdefault('max_n_open_files',
            config['General'].get('max_n_open_files', 300))

    cfg_unzip = config['Unzip']
    def update_from_legacy(new_key, new_section, legacy_key, default=None):
        config.setdefault(new_section, {})
        if new_key in config[new_section]:
            return
        if legacy_key in cfg_unzip:
            config[new_section][new_key] = cfg_unzip[legacy_key]
        elif default is not None:
            config[new_section][new_key] = default
        #else:
        #    msg = 'Please supply "{}" in [{}]'.format(
        #            new_key, new_section)
        #    raise Exception(msg)
    update_from_legacy('JOB_OPTS', 'job.step.unzip.blasr_aln', 'sge_blasr_aln')
    update_from_legacy('JOB_OPTS', 'job.step.unzip.hasm', 'sge_hasm')
    update_from_legacy('JOB_OPTS', 'job.step.unzip.track_reads', 'sge_track_reads')
    update_from_legacy('njobs', 'job.defaults', 'unzip_concurrent_jobs')
    update_from_legacy('njobs', 'job.step.unzip.blasr_aln', 'unzip_blasr_concurrent_jobs')
    update_from_legacy('njobs', 'job.step.unzip.phasing', 'unzip_phasing_concurrent_jobs')
    update_from_legacy('njobs', 'job.step.unzip.quiver', 'quiver_concurrent_jobs')

    if 'smrt_bin' in config['Unzip']:
        LOG.error('You have set option "smrt_bin={}" in the "Unzip" section. That will be ignored. Simply add to your $PATH.'.format(config.get('Unzip', 'smrt_bin')))

    import pprint
    LOG.info('Using config=\n{}'.format(pprint.pformat(config)))

    io.validate_config(config)
    return config

def validate_input_bam_fofn(config, config_fn):
    assert 'input_bam_fofn' in config['Unzip'], 'You must provide "input_bam_fofn" in the [Unzip] section of "{}".'.format(config_fn)
    try:
        # Validate early.
        list(io.yield_bam_fn(config['Unzip']['input_bam_fofn']))
    except Exception:
        msg = 'Failed to validate input_bam_fofn "{}" (from "input_bam_fofn" setting in [Unzip] section of config-file "{}").'.format(
            config['Unzip']['input_bam_fofn'], config_fn)
        LOG.critical(msg)
        raise

def symlink_if_missing(src, name):
    if not os.path.lexists(name):
        LOG.info(' ln -s {} {}'.format(src, name))
        os.symlink(src, name)
    elif os.path.islink(name):
        rn = os.readlink(name)
        if rn == src:
            return # in case 'src' does not exist yet, but they already match
        if not os.path.samefile(src, rn):
            LOG.warning('"{}" != "{}" in {}'.format(src, rn, os.getcwd()))

def update_falcon_symlinks():
    # We might be able to use an older Falcon run if we create the needed symlinks.
    with io.cd('0-rawreads'):
        symlink_if_missing('.', 'build')
        symlink_if_missing('las-gather', 'las-merge-combine')
    with io.cd('1-preads_ovl'):
        symlink_if_missing('.', 'build')
        symlink_if_missing('las-gather', 'las-merge-combine')
    LOG.info('Falcon directories up-to-date.')

def run(target, config_fn, logging_config_fn):
    global LOG
    LOG = support.setup_logger(logging_config_fn)

    config = parse_config(config_fn)
    update_falcon_symlinks()

    # Record the Unzip section as a dict for use by various tasks.
    unzip_config = config['Unzip']
    unzip_config_fn = os.path.join(os.path.dirname(config_fn), 'Unzip_config.json')
    io.serialize(unzip_config_fn, unzip_config) # Some tasks use this.

    if target == 'ccs':
        ccs_all(config, unzip_config_fn)
    else:
        validate_input_bam_fofn(config, config_fn)
        unzip_all(config, unzip_config_fn)
