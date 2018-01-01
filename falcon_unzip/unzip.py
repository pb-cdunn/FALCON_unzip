from falcon_kit import run_support as support
from pypeflow.simple_pwatcher_bridge import (
    PypeLocalFile, makePypeLocalFile, fn,
    PypeTask,
    PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase)
from falcon_kit.FastaReader import FastaReader
from .tasks import unzip as tasks_unzip
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
    unzip_phasing_concurrent_jobs = config['unzip_phasing_concurrent_jobs']
    wf = PypeProcWatcherWorkflow(
        max_jobs=unzip_blasr_concurrent_jobs,
        job_type=config['job_type'],
        job_queue=config.get('job_queue'),
        sge_option=config.get('sge_option'),
        watcher_type=config.get('pwatcher_type'),
        #watcher_directory=config.get('pwatcher_directory', 'mypwatcher'),
        use_tmpdir=config.get('use_tmpdir'),
    )

    read_to_contig_map_file = makePypeLocalFile('3-unzip/reads/get_read_ctg_map/read_to_contig_map')
    # This has lots of inputs from falcon stages 0, 1, and 2.
    wf.addTasks(tasks_unzip.create_tasks_read_to_contig_map(read_to_contig_map_file))

    ctg_list_file = makePypeLocalFile('./3-unzip/reads/ctg_list')
    fofn_file = makePypeLocalFile(config.get('input_fofn', './input.fofn')) # from user config, usually

    track_reads_task = tasks_unzip.get_track_reads_task(config, fofn_file, read_to_contig_map_file, ctg_list_file)
    wf.addTask(track_reads_task)

    # Refresh so that ctg_list_file is available. TODO: Proper scattering.
    wf.refreshTargets()

    gathered_rid_to_phase_file = makePypeLocalFile('./3-unzip/1-hasm/gathered-rid-to-phase/rid_to_phase.all')
    phasing_tasks = list(tasks_unzip.create_phasing_tasks(config, ctg_list_file, gathered_rid_to_phase_file))
    wf.addTasks(phasing_tasks)

    las_fofn_file = makePypeLocalFile('./1-preads_ovl/merge-gather/las.fofn') #'2-asm-falcon/las.fofn'
    job_done = makePypeLocalFile('./3-unzip/1-hasm/hasm_done')
    hasm_task = tasks_unzip.get_hasm_task(config, gathered_rid_to_phase_file, las_fofn_file, job_done)
    wf.addTask(hasm_task)

    wf.max_jobs = unzip_phasing_concurrent_jobs
    wf.refreshTargets()


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
