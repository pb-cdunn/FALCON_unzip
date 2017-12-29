from falcon_kit import run_support as support
from pypeflow.simple_pwatcher_bridge import (
    PypeLocalFile, makePypeLocalFile, fn,
    PypeTask,
    PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase)
from falcon_kit.FastaReader import FastaReader
from .tasks import unzip as tasks_unzip
from .tasks import snakemake
from . import io
import glob
import logging
import os
import re
import time
import ConfigParser

LOG = logging.getLogger(__name__)


def unzip_all(config):
    unzip_blasr_concurrent_jobs = config['unzip_blasr_concurrent_jobs']
    wf = PypeProcWatcherWorkflow(
        max_jobs=unzip_blasr_concurrent_jobs,
        job_type=config['job_type'],
        job_queue=config.get('job_queue'),
        sge_option=config.get('sge_option'),
        watcher_type=config.get('pwatcher_type'),
        #watcher_directory=config.get('pwatcher_directory', 'mypwatcher'),
        use_tmpdir=config.get('use_tmpdir'),
    )
    with open('foo.snake', 'w') as snakemake_writer:
        rule_writer = snakemake.SnakemakeRuleWriter(snakemake_writer)
        tasks_unzip.run_workflow(wf, config, rule_writer)


def run(config_fn, logging_config_fn):
    global LOG
    LOG = support.setup_logger(logging_config_fn)


    config = ConfigParser.ConfigParser()
    config.read(config_fn)

    job_type = 'SGE'
    if config.has_option('General', 'job_type'):
        job_type = config.get('General', 'job_type')

    job_queue = 'default'
    if config.has_option('General', 'job_queue'):
        job_queue = config.get('General', 'job_queue')

    pwatcher_type = 'fs_based'
    if config.has_option('General', 'pwatcher_type'):
        pwatcher_type = config.get('General', 'pwatcher_type')

    sge_blasr_aln = ' -pe smp 24 -q bigmem '
    if config.has_option('Unzip', 'sge_blasr_aln'):
        sge_blasr_aln = config.get('Unzip', 'sge_blasr_aln')

    smrt_bin = ''
    if config.has_option('Unzip', 'smrt_bin'):
        smrt_bin = config.get('Unzip', 'smrt_bin')

    sge_phasing = ' -pe smp 12 -q bigmem'
    if config.has_option('Unzip', 'sge_phasing'):
        sge_phasing = config.get('Unzip', 'sge_phasing')

    sge_hasm = ' -pe smp 48 -q bigmem'
    if config.has_option('Unzip', 'sge_hasm'):
        sge_hasm = config.get('Unzip', 'sge_hasm')

    sge_track_reads = ' -pe smp 12 -q bigmem'
    if config.has_option('Unzip', 'sge_track_reads'):
        sge_track_reads = config.get('Unzip', 'sge_track_reads')

    unzip_blasr_concurrent_jobs = 8
    if config.has_option('Unzip', 'unzip_blasr_concurrent_jobs'):
        unzip_blasr_concurrent_jobs = config.getint('Unzip', 'unzip_blasr_concurrent_jobs')

    unzip_phasing_concurrent_jobs = 8
    if config.has_option('Unzip', 'unzip_phasing_concurrent_jobs'):
        unzip_phasing_concurrent_jobs = config.getint('Unzip', 'unzip_phasing_concurrent_jobs')

    config = {'job_type': job_type,
              'job_queue': job_queue,
              'sge_blasr_aln': sge_blasr_aln,
              'smrt_bin': smrt_bin,
              'sge_phasing': sge_phasing,
              'sge_hasm': sge_hasm,
              'sge_track_reads': sge_track_reads,
              'unzip_blasr_concurrent_jobs': unzip_blasr_concurrent_jobs,
              'unzip_phasing_concurrent_jobs': unzip_phasing_concurrent_jobs,
              'pwatcher_type': pwatcher_type,
              }
    io.update_env_from_config(config, config_fn)

    # support.job_type = 'SGE' #tmp hack until we have a configuration parser

    unzip_all(config)
