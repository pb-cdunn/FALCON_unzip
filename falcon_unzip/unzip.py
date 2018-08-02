from falcon_kit import run_support as support
from falcon_kit import snakemake
from falcon_kit.FastaReader import FastaReader
from pypeflow.simple_pwatcher_bridge import (
    PypeProcWatcherWorkflow,
)
from .tasks import unzip as tasks_unzip
from . import io
import glob
import logging
import os
import re
import time
import falcon_kit.run_support

LOG = logging.getLogger(__name__)


def unzip_all(config):
    job_defaults = config['job.defaults']
    use_tmpdir = job_defaults['use_tmpdir'] # None/False is fine.
    wf = PypeProcWatcherWorkflow(
        job_defaults=job_defaults,
        use_tmpdir=use_tmpdir,
    )
    with open('foo.snake', 'w') as snakemake_writer:
        rule_writer = snakemake.SnakemakeRuleWriter(snakemake_writer)
        tasks_unzip.run_workflow(wf, config, rule_writer)

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
    return config

def parse_config(config_fn):
    config = parse_cfg_file(config_fn)

    # We might need these in the top-level.
    config.setdefault('max_n_open_files',
            config['General'].get('max_n_open_files', 300))

    assert 'input_bam_fofn' in config['Unzip'], 'You must provide "input_bam_fofn" in the [Unzip] section of "{}".'.format(config_fn)

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

def symlink_if_missing(src, name):
    if not os.path.exists(name):
        LOG.info(' ln -s {} {}'.format(src, name))
        os.symlink(src, name)

def update_falcon_symlinks():
    # We might be able to use an older Falcon run if we create the needed symlinks.
    with io.cd('0-rawreads'):
        symlink_if_missing('build', '.')
        symlink_if_missing('las-gather', 'las-merge-combine')
    with io.cd('1-preads_ovl'):
        symlink_if_missing('build', '.')
        symlink_if_missing('las-gather', 'las-merge-combine')
    LOG.info('Falcon directories up-to-date.')

def run(config_fn, logging_config_fn):
    global LOG
    LOG = support.setup_logger(logging_config_fn)

    config = parse_config(config_fn)
    update_falcon_symlinks()

    unzip_all(config)
