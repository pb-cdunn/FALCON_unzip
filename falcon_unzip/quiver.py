from falcon_kit import run_support as support
#from falcon_kit import snakemake
from .tasks import snakemake
from pypeflow.simple_pwatcher_bridge import (
    PypeProcWatcherWorkflow)
from .tasks import quiver as tasks_quiver
from . import io
import glob
import logging
import os
import re
import time
import ConfigParser

LOG = logging.getLogger(__name__)


def run(config_fn, logging_config_fn):
    global LOG
    LOG = support.setup_logger(logging_config_fn)

    config_absbasedir = os.path.dirname(os.path.abspath(config_fn))

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

    max_n_open_files = 300
    if config.has_option('General', 'max_n_open_files'):
        max_n_open_files = config.getint('General', 'max_n_open_files')

    sge_option = ''
    if config.has_option('Unzip', 'sge_option'):
        sge_option = config.get('Unzip', 'sge_option')

    sge_track_reads = sge_option
    if config.has_option('Unzip', 'sge_track_reads'):
        sge_track_reads = config.get('Unzip', 'sge_track_reads')

    sge_quiver = sge_option
    if config.has_option('Unzip', 'sge_quiver'):
        sge_quiver = config.get('Unzip', 'sge_quiver')

    input_bam_fofn = 'input_bam.fofn'
    if config.has_option('Unzip', 'input_bam_fofn'):
        input_bam_fofn = config.get('Unzip', 'input_bam_fofn')
    if not os.path.isabs(input_bam_fofn):
        input_bam_fofn = os.path.join(config_absbasedir, input_bam_fofn)

    quiver_concurrent_jobs = 8
    if config.has_option('Unzip', 'quiver_concurrent_jobs'):
        quiver_concurrent_jobs = config.getint('Unzip', 'quiver_concurrent_jobs')

    if config.has_option('Unzip', 'smrt_bin'):
        LOG.warning('You have set option "smrt_bin={}" in the "Unzip" section. That will be ignored. Simply and to your $PATH.'.format(config.get('Unzip', 'smrt_bin')))

    config = {'job_type': job_type,
              'job_queue': job_queue,
              'sge_option': sge_option,
              'sge_quiver': sge_quiver,
              'sge_track_reads': sge_track_reads,
              'input_bam_fofn': input_bam_fofn,
              'max_n_open_files': max_n_open_files,
              'pwatcher_type': pwatcher_type,
    }
    io.update_env_from_config(config, config_fn)

    # support.job_type = 'SGE' #tmp hack until we have a configuration parser

    wf = PypeProcWatcherWorkflow(
        max_jobs=quiver_concurrent_jobs,
        job_type=config['job_type'],
        job_queue=config.get('job_queue'),
        sge_option=config.get('sge_option'),
        watcher_type=config.get('pwatcher_type'),
        #watcher_directory=config.get('pwatcher_directory', 'mypwatcher'),
        use_tmpdir=config.get('use_tmpdir'),
    )
    with open('foo.snake', 'w') as snakemake_writer:
        rule_writer = snakemake.SnakemakeRuleWriter(snakemake_writer)
        tasks_quiver.run_workflow(wf, config, rule_writer)
