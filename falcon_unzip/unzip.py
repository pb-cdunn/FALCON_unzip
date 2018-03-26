from falcon_kit import run_support as support
from falcon_kit import snakemake
#from .tasks import snakemake
from falcon_kit.FastaReader import FastaReader
from pypeflow.simple_pwatcher_bridge import (
    PypeProcWatcherWorkflow)
from .tasks import unzip as tasks_unzip
from . import io
import glob
import logging
import os
import re
import time
#import ConfigParser
import falcon_kit.run_support

LOG = logging.getLogger(__name__)


def unzip_all(config):
    wf = PypeProcWatcherWorkflow(
        job_defaults=config['job.defaults'],
    )
    with open('foo.snake', 'w') as snakemake_writer:
        rule_writer = snakemake.SnakemakeRuleWriter(snakemake_writer)
        tasks_unzip.run_workflow(wf, config, rule_writer)

def update_config_from_sections(config, cfg):
    allowed_sections = set([
            'General',
            'Unzip',
            'job.step.uzip.track_reads',
            'job.step.unzip.phasing',
            'job.step.unzip.blasr_aln',
            'job.step.unzip.hasm',
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
    """
    # New: Parse sections (case-sensitively), into sub-dicts.
    config = dict()
    with open(config_fn) as stream:
        cfg2 = falcon_kit.run_support.parse_cfg_with_sections(stream)
        update_config_from_sections(config, cfg2)
    falcon_kit.run_support.update_job_defaults_section(config)
    return config

def run(config_fn, logging_config_fn):
    global LOG
    LOG = support.setup_logger(logging_config_fn)

    config = parse_cfg_file(config_fn)

    cfg_unzip = config['Unzip']

    def update_from_legacy(new_key, new_section, legacy_key, default=None):
        config.setdefault(new_section, {})
        if legacy_key in cfg_unzip and new_key not in config[new_section]:
            config[new_section][new_key] = cfg_unzip[legacy_key]
        elif default is not None:
            config[new_section][new_key] = default
    update_from_legacy('JOB_OPTS', 'job.step.unzip.blasr_aln', 'sge_blasr_aln')
    update_from_legacy('JOB_OPTS', 'job.step.unzip.hasm', 'sge_hasm')
    update_from_legacy('JOB_OPTS', 'job.step.unzip.track_reads', 'sge_track_reads')
    update_from_legacy('njobs', 'job.defaults', 'unzip_blasr_concurrent_jobs', default=8)
    update_from_legacy('njobs', 'job.step.unzip.blasr_aln', 'unzip_blasr_concurrent_jobs', default=8)
    update_from_legacy('njobs', 'job.step.unzip.phasing', 'unzip_phasing_concurrent_jobs', default=8)

    if 'smrt_bin' in config['Unzip']:
        LOG.error('You have set option "smrt_bin={}" in the "Unzip" section. That will be ignored. Simply add to your $PATH.'.format(config.get('Unzip', 'smrt_bin')))

    unzip_all(config)
