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

def task_track_reads(self):
    job_done = fn(self.job_done)
    fofn_fn = os.path.relpath(fn(self.fofn))

    topdir = os.path.relpath(self.parameters['topdir'])

    script_fn = 'track_reads.sh'
    script = """\
set -vex
trap 'touch {job_done}.exit' EXIT
hostname
date

python -m falcon_unzip.mains.rr_ctg_track --base-dir={topdir} --output=rawread_to_contigs
python -m falcon_unzip.mains.pr_ctg_track --base-dir={topdir} --output=pread_to_contigs
# Those outputs are used only by fetch_reads.
python -m falcon_unzip.mains.fetch_reads --base-dir={topdir} --fofn={fofn_fn}
date
touch {job_done}
""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def task_hasm(self):
    rid_to_phase_all = fn(self.rid_to_phase_all)
    las_fofn = fn(self.las_fofn)
    job_done = fn(self.job_done)
    #config = self.parameters['config']

    script_fn = 'hasm.sh'
    script = """\
set -vex
trap 'touch {job_done}.exit' EXIT
hostname
date

python -m falcon_unzip.mains.ovlp_filter_with_phase --fofn {las_fofn} --max_diff 120 --max_cov 120 --min_cov 1 --n_core 48 --min_len 2500 --db ../../1-preads_ovl/preads.db --rid_phase_map {rid_to_phase_all} > preads.p_ovl
python -m falcon_unzip.mains.phased_ovlp_to_graph preads.p_ovl --min_len 2500 > fc.log
if [ -e ../../1-preads_ovl/preads4falcon.fasta ];
then
  ln -sf ../../1-preads_ovl/preads4falcon.fasta .
else
  ln -sf ../../1-preads_ovl/db2falcon/preads4falcon.fasta .
fi
python -m falcon_unzip.mains.graphs_to_h_tigs --fc_asm_path ../../2-asm-falcon/ --fc_hasm_path ./ --ctg_id all --rid_phase_map {rid_to_phase_all} --fasta preads4falcon.fasta

# more script -- a little bit hacky here, we should improve

WD=$PWD
for f in `cat ../reads/ctg_list `; do mkdir -p $WD/$f; cd $WD/$f; python -m falcon_unzip.mains.dedup_h_tigs $f; done

## prepare for quviering the haplotig
cd $WD/..

find 0-phasing -name "phased_reads" | sort | xargs cat >| all_phased_reads
find 1-hasm -name "h_ctg_ids.*" | sort | xargs cat >| all_h_ctg_ids
find 1-hasm -name "p_ctg_edges.*" | sort | xargs cat >| all_p_ctg_edges
find 1-hasm -name "h_ctg_edges.*" | sort | xargs cat >| all_h_ctg_edges
find 1-hasm -name "p_ctg.*.fa" | sort | xargs cat >| all_p_ctg.fa
find 1-hasm -name "h_ctg.*.fa" | sort | xargs cat >| all_h_ctg.fa

# Generate a GFA for only primary contigs and haplotigs.
time python -m falcon_unzip.mains.unzip_gen_gfa_v1 --unzip-root $WD/.. --p-ctg-fasta $WD/../all_p_ctg.fa --h-ctg-fasta $WD/../all_h_ctg.fa --preads-fasta $WD/preads4falcon.fasta >| $WD/../asm.gfa

# Generate a GFA of all assembly graph edges. This GFA can contain
# edges and nodes which are not part of primary contigs and haplotigs
time python -m falcon_unzip.mains.unzip_gen_gfa_v1 --unzip-root $WD/.. --p-ctg-fasta $WD/../all_p_ctg.fa --h-ctg-fasta $WD/../all_h_ctg.fa --preads-fasta $WD/preads4falcon.fasta --add-string-graph >| $WD/../sg.gfa

cd ../
date
touch {job_done}
""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


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
    wf.addTasks(tasks_unzip.create_tasks_read_to_contig_map(read_to_contig_map_file))

    ctg_list_file = makePypeLocalFile('./3-unzip/reads/ctg_list')
    fofn_file = makePypeLocalFile('./input.fofn') # TODO: Make explicit input from user.

    parameters = {'config': config,
                  'sge_option': config['sge_track_reads'],
                  'topdir': os.getcwd(),
                  }
    job_done = makePypeLocalFile('./3-unzip/reads/track_reads_done')
    make_track_reads_task = PypeTask(
            inputs={
                'fofn': fofn_file,
                'read_to_contig_map': read_to_contig_map_file,
            },
            outputs={
                'job_done': job_done, 'ctg_list_file': ctg_list_file,
            },
            parameters=parameters,
            )
    track_reads_task = make_track_reads_task(task_track_reads)
    wf.addTask(track_reads_task)

    # For refresh so that ctg_list_file is available. TODO: Proper scattering.
    wf.refreshTargets()

    gathered_rid_to_phase_file = makePypeLocalFile('./3-unzip/1-hasm/gathered-rid-to-phase/rid_to_phase.all')
    phasing_tasks = list(tasks_unzip.create_phasing_tasks(config, ctg_list_file, gathered_rid_to_phase_file))
    wf.addTasks(phasing_tasks)

    parameters['sge_option'] = config['sge_hasm']
    job_done = makePypeLocalFile('./3-unzip/1-hasm/hasm_done')
    make_hasm_task = PypeTask(
            inputs={
                'rid_to_phase_all': gathered_rid_to_phase_file,
                'las_fofn': makePypeLocalFile('./1-preads_ovl/merge-gather/las.fofn'), #'2-asm-falcon/las.fofn',
            },
            outputs={
                'job_done': job_done,
            },
            parameters=parameters,
            )
    hasm_task = make_hasm_task(task_hasm)
    wf.addTask(hasm_task)

    wf.max_jobs = unzip_phasing_concurrent_jobs
    wf.refreshTargets()


def run(config_fn):
    global LOG
    LOG = support.setup_logger(None)

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
