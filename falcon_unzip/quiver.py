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
    }
    input_bam_fofn_fn = config['input_bam_fofn']
    input_bam_fofn_plf = makePypeLocalFile(input_bam_fofn_fn)
    hasm_done_plf = makePypeLocalFile('./3-unzip/1-hasm/hasm_done')  # by convention

    track_reads_h_done_plf = makePypeLocalFile('./4-quiver/track_reads/track_reads_h_done')
    make_task = PypeTask(inputs={
        'input_bam_fofn': input_bam_fofn_plf,
        'hasm_done': hasm_done_plf},
        outputs={'job_done': track_reads_h_done_plf},
        parameters=parameters,
    )
    wf.addTask(make_task(tasks_quiver.task_track_reads_h))
    # Note: The output is actually './4-quiver/track_reads/rawread_to_contigs

    read2ctg_plf = makePypeLocalFile('./4-quiver/select_reads/read2ctg.msgpack')
    make_task = PypeTask(inputs={
        # Some implicit inputs, plus these deps:
        'track_reads_h_done': track_reads_h_done_plf,
        'input_bam_fofn': input_bam_fofn_plf,
        'hasm_done': hasm_done_plf},
        outputs={
        'read2ctg': read2ctg_plf},
        parameters=parameters,
    )
    wf.addTask(make_task(tasks_quiver.task_select_reads_h))

    merged_fofn_plf = makePypeLocalFile('./4-quiver/merge_reads/merged.fofn')
    make_task = PypeTask(inputs={
        'input_bam_fofn': input_bam_fofn_plf,
        'read2ctg': read2ctg_plf},
        outputs={
        'merged_fofn': merged_fofn_plf},
        parameters=parameters,
    )
    wf.addTask(make_task(tasks_quiver.task_merge_reads))

    scattered_segregate_plf = makePypeLocalFile('./4-quiver/segregate_scatter/scattered.json')
    make_task = PypeTask(
        inputs={
            'merged_fofn': merged_fofn_plf,
        },
        outputs={
            'scattered_segregate_json': scattered_segregate_plf,
        },
        parameters=parameters,
    )
    wf.addTask(make_task(tasks_quiver.task_segregate_scatter))
    wf.refreshTargets()

    # Segregate reads from merged BAM files in parallel.
    # (If this were not done in Python, it could probably be in serial.)
    jn2segregated_bam_fofn = tasks_quiver.create_segregate_jobs(wf, parameters, scattered_segregate_plf)
    # ctg is encoded into each filepath within each FOFN.

    ctg2segregated_bamfn_plf = makePypeLocalFile('./4-quiver/segregate_gather/ctg2segregated_bamfn.msgpack')
    make_task = PypeTask(
        inputs=jn2segregated_bam_fofn,
        outputs={
            'ctg2segregated_bamfn': ctg2segregated_bamfn_plf,
        },
        parameters=parameters,
    )
    wf.addTask(make_task(tasks_quiver.task_segregate_gather))
    wf.refreshTargets()

    scattered_quiver_plf = makePypeLocalFile('4-quiver/quiver_scatter/scattered.json')
    parameters = {
        'config': config,
    }
    make_task = PypeTask(
        inputs={
            'p_ctg_fa': makePypeLocalFile('3-unzip/all_p_ctg.fa'),
            'h_ctg_fa': makePypeLocalFile('3-unzip/all_h_ctg.fa'),
            'ctg2bamfn': ctg2segregated_bamfn_plf,
        },
        outputs={
            'scattered_quiver_json': scattered_quiver_plf,
        },
        parameters=parameters,
    )
    wf.addTask(make_task(tasks_quiver.task_scatter_quiver))
    wf.refreshTargets()

    p_ctg_out, h_ctg_out, job_done_plfs = tasks_quiver.create_quiver_jobs(wf, scattered_quiver_plf)

    gathered_p_ctg_plf = makePypeLocalFile('4-quiver/cns_gather/p_ctg.txt')
    gathered_h_ctg_plf = makePypeLocalFile('4-quiver/cns_gather/h_ctg.txt')
    gather_done_plf = makePypeLocalFile('4-quiver/cns_gather/job_done')
    io.mkdirs('4-quiver/cns_gather')
    with open(fn(gathered_p_ctg_plf), 'w') as ifs:
        for cns_fasta_fn, cns_fastq_fn in sorted(p_ctg_out):
            ifs.write('{} {}\n'.format(cns_fasta_fn, cns_fastq_fn))
    with open(fn(gathered_h_ctg_plf), 'w') as ifs:
        for cns_fasta_fn, cns_fastq_fn in sorted(h_ctg_out):
            ifs.write('{} {}\n'.format(cns_fasta_fn, cns_fastq_fn))

    make_task = PypeTask(
        inputs=job_done_plfs,
        outputs={
            'job_done': gather_done_plf,
        },
        parameters={},
    )
    wf.addTask(make_task(tasks_quiver.task_gather_quiver))
    wf.refreshTargets()

    cns_p_ctg_fasta_plf = makePypeLocalFile('4-quiver/cns_output/cns_p_ctg.fasta')
    cns_p_ctg_fastq_plf = makePypeLocalFile('4-quiver/cns_output/cns_p_ctg.fastq')
    cns_h_ctg_fasta_plf = makePypeLocalFile('4-quiver/cns_output/cns_h_ctg.fasta')
    cns_h_ctg_fastq_plf = makePypeLocalFile('4-quiver/cns_output/cns_h_ctg.fastq')
    zcat_done_plf = makePypeLocalFile('4-quiver/cns_output/job_done')
    make_task = PypeTask(
        inputs={
            'gathered_p_ctg': gathered_p_ctg_plf,
            'gathered_h_ctg': gathered_h_ctg_plf,
            'gather_done': gather_done_plf,
        },
        outputs={
            'cns_p_ctg_fasta': cns_p_ctg_fasta_plf,
            'cns_p_ctg_fastq': cns_p_ctg_fastq_plf,
            'cns_h_ctg_fasta': cns_h_ctg_fasta_plf,
            'cns_h_ctg_fastq': cns_h_ctg_fastq_plf,
            'job_done': zcat_done_plf,
        },
    )
    wf.addTask(make_task(tasks_quiver.task_cns_zcat))

    wf.refreshTargets()
