from falcon_kit import run_support as support
from pypeflow.simple_pwatcher_bridge import (
    makePypeLocalFile, fn,
    PypeTask,
    PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase)
from .tasks import quiver as tasks_quiver
from . import io
import glob
import logging
import os
import re
import time
import ConfigParser

LOG = logging.getLogger(__name__)


def run(config_fn):
    global LOG
    LOG = support.setup_logger(None)

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

    sge_track_reads = ' -pe smp 12 -q bigmem'
    if config.has_option('Unzip', 'sge_track_reads'):
        sge_track_reads = config.get('Unzip', 'sge_track_reads')

    sge_quiver = ' -pe smp 24 -q bigmem '
    if config.has_option('Unzip', 'sge_quiver'):
        sge_quiver = config.get('Unzip', 'sge_quiver')

    smrt_bin = ''
    if config.has_option('Unzip', 'smrt_bin'):
        smrt_bin = config.get('Unzip', 'smrt_bin')

    input_bam_fofn = 'input_bam.fofn'
    if config.has_option('Unzip', 'input_bam_fofn'):
        input_bam_fofn = config.get('Unzip', 'input_bam_fofn')
    if not os.path.isabs(input_bam_fofn):
        input_bam_fofn = os.path.join(config_absbasedir, input_bam_fofn)

    quiver_concurrent_jobs = 8
    if config.has_option('Unzip', 'quiver_concurrent_jobs'):
        quiver_concurrent_jobs = config.getint('Unzip', 'quiver_concurrent_jobs')

    config = {'job_type': job_type,
              'job_queue': job_queue,
              'sge_quiver': sge_quiver,
              'sge_track_reads': sge_track_reads,
              'input_bam_fofn': input_bam_fofn,
              'max_n_open_files': max_n_open_files,
              'pwatcher_type': pwatcher_type,
              'smrt_bin': smrt_bin}
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

    abscwd = os.path.abspath('.')
    parameters = {
        'sge_option': config['sge_track_reads'],  # applies to select_reads task also, for now
        'max_n_open_files': config['max_n_open_files'],
        'topdir': os.getcwd(),
    }
    input_bam_fofn_fn = config['input_bam_fofn']
    input_bam_fofn_plf = makePypeLocalFile(input_bam_fofn_fn)
    hasm_done_plf = makePypeLocalFile('./3-unzip/1-hasm/hasm_done')  # by convention
    track_reads_h_done_plf = makePypeLocalFile('./4-quiver/track_reads/track_reads_h_done')
    track_reads_rr2c_plf = makePypeLocalFile('./4-quiver/track_reads/rawread_to_contigs')
    wf.addTask(tasks_quiver.get_track_reads_h_task(
        parameters, input_bam_fofn_plf, hasm_done_plf,
        track_reads_h_done_plf, track_reads_rr2c_plf))

    read2ctg_plf = makePypeLocalFile('./4-quiver/select_reads/read2ctg.msgpack')
    wf.addTask(tasks_quiver.get_select_reads_h_task(
        parameters, track_reads_h_done_plf, input_bam_fofn_plf, hasm_done_plf,
        read2ctg_plf))

    merged_fofn_plf = makePypeLocalFile('./4-quiver/merge_reads/merged.fofn')
    wf.addTask(tasks_quiver.get_merge_reads_task(
        parameters, input_bam_fofn_plf, read2ctg_plf, merged_fofn_plf))

    scattered_segregate_plf = makePypeLocalFile('./4-quiver/segregate_scatter/scattered.json')
    wf.addTask(tasks_quiver.get_segregate_scatter_task(
        parameters, merged_fofn_plf, scattered_segregate_plf))
    wf.refreshTargets()

    ctg2segregated_bamfn_plf = makePypeLocalFile('./4-quiver/segregate_gather/ctg2segregated_bamfn.msgpack')
    wf.addTasks(list(tasks_quiver.yield_segregate_bam_tasks(
        parameters, scattered_segregate_plf, ctg2segregated_bamfn_plf)))

    scattered_quiver_plf = makePypeLocalFile('4-quiver/quiver_scatter/scattered.json')
    parameters = {
        'config': config,
    }
    wf.addTask(tasks_quiver.get_scatter_quiver_task(
        parameters, ctg2segregated_bamfn_plf,
        scattered_quiver_plf,
        ))
    wf.refreshTargets()

    gathered_p_ctg_plf = makePypeLocalFile('4-quiver/cns_gather/p_ctg.txt')
    gathered_h_ctg_plf = makePypeLocalFile('4-quiver/cns_gather/h_ctg.txt')

    wf.addTasks(list(tasks_quiver.yield_quiver_tasks(
        scattered_quiver_plf,
        gathered_p_ctg_plf, gathered_h_ctg_plf)))

    zcat_done_plf = makePypeLocalFile('4-quiver/cns_output/job_done')

    wf.addTask(tasks_quiver.get_cns_zcat_task(
        gathered_p_ctg_plf, gathered_h_ctg_plf, gather_done_plf,
        zcat_done_plf))

    wf.refreshTargets()
