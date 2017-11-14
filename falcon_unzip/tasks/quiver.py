from pypeflow.simple_pwatcher_bridge import (
    makePypeLocalFile, fn,
    PypeTask,
)
from pypeflow.sample_tasks import gen_task # TODO: Better path
from falcon_kit.FastaReader import FastaReader
from .. import io
import json
import logging
import os
import re

LOG = logging.getLogger(__name__)

# For now, in/outputs are in various directories, by convention, including '0-rawreads/m_*/*.msgpack'
TASK_TRACK_READS_H_SCRIPT = """\
python -m falcon_unzip.mains.get_read_hctg_map --base-dir={topdir} --output=read_to_contig_map
# formerly generated ./4-quiver/read_maps/read_to_contig_map

fc_rr_hctg_track.py --base-dir={topdir} --stream
# That writes into 0-rawreads/m_*/

abs_rawread_to_contigs=$(readlink -f {rawread_to_contigs}) #TODO: No readlink
cd {topdir}
fc_rr_hctg_track2.exe --output=${{abs_rawread_to_contigs}}
cd -
ls -l {rawread_to_contigs}
"""


# For now, in/outputs are in various directories, by convention.
TASK_SELECT_READS_H_SCRIPT = """\
python -m falcon_unzip.mains.get_read2ctg --output={read2ctg} {input_bam_fofn}
"""


def task_merge_reads(self):
    merged_fofn = fn(self.merged_fofn)
    read2ctg = fn(self.read2ctg)
    input_bam_fofn = fn(self.input_bam_fofn)
    max_n_open_files = self.parameters['max_n_open_files']
    topdir = os.path.relpath(self.parameters['topdir'])
    script_fn = 'merge_reads.sh'

    # For now, in/outputs are in various directories, by convention.
    script = """\
set -vex
#trap 'touch {merged_fofn}.exit' EXIT
hostname
date

cd {topdir}
#fc_select_reads_from_bam.py --max-n-open-files={max_n_open_files} {input_bam_fofn}
pwd
python -m falcon_unzip.mains.bam_partition_and_merge --max-n-open-files={max_n_open_files} --read2ctg-fn={read2ctg} --merged-fn={merged_fofn} {input_bam_fofn}

date
cd -
# Expect {merged_fofn}
""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def task_run_quiver(self):
    ref_fasta = fn(self.ref_fasta)
    read_bam = fn(self.read_bam)

    cns_fasta = fn(self.cns_fasta)
    cns_fastq = fn(self.cns_fastq)
    job_done = fn(self.job_done)

    job_uid = self.parameters['job_uid']
    ctg_id = self.parameters['ctg_id']

    # TODO: tmpdir

    script_fn = 'cns_%s.sh' % (ctg_id)
    script = """\
set -vex
trap 'touch {job_done}.exit' EXIT
hostname
date

samtools faidx {ref_fasta}
pbalign --tmpDir=/localdisk/scratch/ --nproc=24 --minAccuracy=0.75 --minLength=50\
          --minAnchorSize=12 --maxDivergence=30 --concordant --algorithm=blasr\
          --algorithmOptions=--useQuality --maxHits=1 --hitPolicy=random --seed=1\
            {read_bam} {ref_fasta} aln-{ctg_id}.bam
#python -c 'import ConsensusCore2 as cc2; print cc2' # So quiver likely works.
(variantCaller --algorithm=arrow -x 5 -X 120 -q 20 -j 24 -r {ref_fasta} aln-{ctg_id}.bam\
            -o {cns_fasta} -o {cns_fastq}) || echo WARNING quiver failed. Maybe no reads for this block.
date
touch {job_done}
""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def task_cns_zcat(self):
    gathered_quiver = fn(self.gathered_quiver)
    cns_p_ctg_fasta = fn(self.cns_p_ctg_fasta)
    cns_p_ctg_fastq = fn(self.cns_p_ctg_fastq)
    cns_h_ctg_fasta = fn(self.cns_h_ctg_fasta)
    cns_h_ctg_fastq = fn(self.cns_h_ctg_fastq)

    script_fn = 'cns_zcat.sh'
    script = """\
python -m falcon_unzip.mains.cns_zcat \
    --gathered-quiver-fn={gathered_quiver} \
    --cns-p-ctg-fasta-fn={cns_p_ctg_fasta} \
    --cns-p-ctg-fastq-fn={cns_p_ctg_fastq} \
    --cns-h-ctg-fasta-fn={cns_h_ctg_fasta} \
    --cns-h-ctg-fastq-fn={cns_h_ctg_fastq} \

""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def task_scatter_quiver(self):
    p_ctg_fn = fn(self.p_ctg_fa)
    h_ctg_fn = fn(self.h_ctg_fa)
    out_json = fn(self.scattered_quiver_json)
    ctg2bamfn_fn = fn(self.ctg2bamfn)
    config = self.parameters['config']

    ctg2bamfn = io.deserialize(ctg2bamfn_fn)

    ref_seq_data = {}

    # I think this will crash if the file is empty. Maybe that is ok.
    p_ctg_fa = FastaReader(p_ctg_fn)
    ctg_types = {}
    for r in p_ctg_fa:
        rid = r.name.split()[0]
        ref_seq_data[rid] = r.sequence
        ctg_types[rid] = 'p'

    # I think this will crash if the file is empty. Maybe that is ok.
    h_ctg_fa = FastaReader(h_ctg_fn)
    for r in h_ctg_fa:
        rid = r.name.split()[0]
        ref_seq_data[rid] = r.sequence
        ctg_types[rid] = 'h'

    ctg_ids = sorted(ref_seq_data.keys())
    # p_ctg_out=[]
    # h_ctg_out=[]
    #job_done_plfs = {}
    jobs = []
    for ctg_id in ctg_ids:
        sequence = ref_seq_data[ctg_id]
        m_ctg_id = ctg_id.split('-')[0]
        wd = os.path.join(os.getcwd(), m_ctg_id)
        ref_fasta = os.path.join(wd, '{ctg_id}_ref.fa'.format(ctg_id=ctg_id))
        #cns_fasta = makePypeLocalFile(os.path.join(wd, 'cns-{ctg_id}.fasta.gz'.format(ctg_id = ctg_id)))
        #cns_fastq = makePypeLocalFile(os.path.join(wd, 'cns-{ctg_id}.fastq.gz'.format(ctg_id = ctg_id)))
        #job_done = makePypeLocalFile(os.path.join(wd, '{ctg_id}_quiver_done'.format(ctg_id = ctg_id)))
        ctg_types2 = {}
        ctg_types2[ctg_id] = ctg_types[ctg_id]

        # if os.path.exists(read_bam):
        if ctg_id in ctg2bamfn:
            read_bam = ctg2bamfn[ctg_id]
            # The segregated *.sam are created in task_segregate.
            # Network latency should not matter because we have already waited for the 'done' file.
            io.mkdirs(wd)
            if not os.path.exists(ref_fasta):
                # TODO(CD): Up to 50MB of seq data. Should do this on remote host.
                #   See https://github.com/PacificBiosciences/FALCON_unzip/issues/59
                with open(ref_fasta, 'w') as f:
                    print >>f, '>' + ctg_id
                    print >>f, sequence
            new_job = {}
            new_job['ctg_id'] = ctg_id
            new_job['ctg_types'] = ctg_types2
            new_job['smrt_bin'] = config['smrt_bin']
            new_job['sge_option'] = config['sge_quiver']
            new_job['ref_fasta'] = ref_fasta
            new_job['read_bam'] = read_bam
            jobs.append(new_job)
    io.serialize(out_json, jobs)


def task_gather_quiver(self):
    """We wrote the "gathered" files during task construction.
    """
    p_ctg_out = list()
    h_ctg_out = list()
    re_done = re.compile(r'([^/]*)_([ph])_quiver_done')
    for k,v in self.inputs.iteritems():
        print 'items:', k, v
        """
        m_ctg_id = ctg_id.split('-')[0]
        wd = os.path.join(os.getcwd(), './4-quiver/', m_ctg_id)
        #ref_fasta = makePypeLocalFile(os.path.join(wd, '{ctg_id}_ref.fa'.format(ctg_id = ctg_id)))
        #read_bam = makePypeLocalFile(os.path.join(os.getcwd(), './4-quiver/reads/' '{ctg_id}.sam'.format(ctg_id = ctg_id)))
        cns_fasta = makePypeLocalFile(os.path.join(wd, 'cns-{ctg_id}.fasta.gz'.format(ctg_id=ctg_id)))
        cns_fastq = makePypeLocalFile(os.path.join(wd, 'cns-{ctg_id}.fastq.gz'.format(ctg_id=ctg_id)))
        job_done = makePypeLocalFile(os.path.join(wd, '{ctg_id}_quiver_done'.format(ctg_id=ctg_id)))
        """
        wd = os.path.dirname(fn(v))
        basename = os.path.basename(fn(v))
        mo = re_done.search(basename)
        if not mo:
            raise Exception('No match: {!r} not in {!r}'.format(
                basename, re_done.pattern))
        ctg_id = mo.group(1)
        ctg_type = mo.group(2)
        cns_fasta = '{wd}/cns-{ctg_id}.fasta.gz'.format(**locals())
        cns_fastq = '{wd}/cns-{ctg_id}.fastq.gz'.format(**locals())
        assert ctg_type in 'ph', 'ctg_type={!r}'.format(ctg_type)
        if ctg_type == 'p':
            p_ctg_out.append([cns_fasta, cns_fastq])
        elif ctg_type == 'h':
            h_ctg_out.append([cns_fasta, cns_fastq])
    gathered_quiver = fn(self.gathered_quiver)
    gathered_quiver_dict = {
            'p_ctg': list(sorted(p_ctg_out)), # cns_fasta_fn, cns_fastq_fn
            'h_ctg': list(sorted(h_ctg_out)), # cns_fasta_fn, cns_fastq_fn
    }
    io.serialize(gathered_quiver, gathered_quiver_dict)



def task_segregate_scatter(self):
    merged_fofn_fn = fn(self.merged_fofn)
    scattered_segregate_json_fn = fn(self.scattered_segregate_json)

    LOG.info('Scatting segregate-reads tasks. Reading merged BAM names from FOFN: {!r}'.format(
        merged_fofn_fn))
    fns = list(io.yield_abspath_from_fofn(merged_fofn_fn))
    jobs = list()
    for i, merged_bamfn in enumerate(fns):
        job = dict()
        job['merged_bamfn'] = merged_bamfn
        job_name = 'segr{:03d}'.format(i)
        job['job_name'] = job_name
        jobs.append(job)

    io.serialize(scattered_segregate_json_fn, jobs)
    # Fast (for now), so do it locally.


def task_run_segregate(self):
    # max_n_open_files = 300 # Ignored for now. Should not matter here.
    merged_bamfn = fn(self.merged_bamfn)
    segregated_bam_fofn = fn(self.segregated_bam_fofn)

    script = """
python -m falcon_unzip.mains.bam_segregate --merged-fn={merged_bamfn} --output-fn={segregated_bam_fofn}
""".format(**locals())

    script_fn = 'run_bam_segregate.sh'
    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def task_segregate_gather(self):
    jn2segregated_bam_fofn = self.inputs
    ctg2segregated_bamfn_fn = fn(self.ctg2segregated_bamfn)

    ctg2segregated_bamfn = dict()
    for jn, plf in jn2segregated_bam_fofn.iteritems():
        # We do not really care about the arbitrary job-name.
        fofn_fn = fn(plf)
        # Read FOFN.
        segregated_bam_fns = list(io.yield_abspath_from_fofn(fofn_fn))
        # Discern ctgs from filepaths.
        for bamfn in segregated_bam_fns:
            basename = os.path.basename(bamfn)
            ctg = os.path.splitext(basename)[0]
            ctg2segregated_bamfn[ctg] = bamfn
    io.serialize(ctg2segregated_bamfn_fn, ctg2segregated_bamfn)
    io.serialize(ctg2segregated_bamfn_fn + '.json', ctg2segregated_bamfn)  # for debugging
    # Do not generate a script. This is light and fast, so do it locally.


def yield_segregate_bam_tasks(parameters, scattered_segregate_plf, ctg2segregated_bamfn_plf):
    # Segregate reads from merged BAM files in parallel.
    # (If this were not done in Python, it could probably be in serial.)

    jn2segregated_bam_fofn = dict()  # job_name -> FOFN_plf
    # ctg is encoded into each filepath within each FOFN.

    scattered_segregate_fn = fn(scattered_segregate_plf)
    jobs = io.deserialize(scattered_segregate_fn)
    basedir = os.path.dirname(scattered_segregate_fn)  # Should this be relative to cwd?
    for job in jobs:
        job_name = job['job_name']
        merged_bamfn_plf = makePypeLocalFile(job['merged_bamfn'])
        wd = os.path.join(basedir, job_name)
        # ctg is encoded into each filepath within the FOFN.
        segregated_bam_fofn_plf = makePypeLocalFile(os.path.join(wd, 'segregated_bam.fofn'))
        make_task = PypeTask(
            inputs={
                # The other input is next to this one, named by convention.
                'merged_bamfn': merged_bamfn_plf},
            outputs={
                'segregated_bam_fofn': segregated_bam_fofn_plf},
            parameters=parameters,
        )
        yield make_task(task_run_segregate)
        jn2segregated_bam_fofn[job_name] = segregated_bam_fofn_plf

    make_task = PypeTask(
        inputs=jn2segregated_bam_fofn,
        outputs={
            'ctg2segregated_bamfn': ctg2segregated_bamfn_plf,
        },
        parameters=parameters,
    )
    yield make_task(task_segregate_gather)


def get_scatter_quiver_task(
        parameters, ctg2segregated_bamfn_plf,
        scattered_quiver_plf,
):
    make_task = PypeTask(
        inputs={
            'p_ctg_fa': makePypeLocalFile('3-unzip/all_p_ctg.fa'), # TODO: make explicit
            'h_ctg_fa': makePypeLocalFile('3-unzip/all_h_ctg.fa'),
            'ctg2bamfn': ctg2segregated_bamfn_plf,
        },
        outputs={
            'scattered_quiver_json': scattered_quiver_plf,
        },
        parameters=parameters,
    )
    return make_task(task_scatter_quiver)

def yield_quiver_tasks(
        scattered_quiver_plf,
        gathered_quiver_plf,
):
    scattered_quiver_fn = fn(scattered_quiver_plf)
    jobs = json.loads(open(scattered_quiver_fn).read())
    #ctg_ids = sorted(jobs['ref_seq_data'])
    p_ctg_out = []
    h_ctg_out = []
    job_done_plfs = {}
    for job in jobs:
        ctg_id = job['ctg_id']
        ctg_types = job['ctg_types']
        smrt_bin = job['smrt_bin']
        sge_option = job['sge_option']
        ref_fasta = makePypeLocalFile(job['ref_fasta'])
        read_bam = makePypeLocalFile(job['read_bam'])
        m_ctg_id = ctg_id.split('-')[0]
        wd = os.path.join(os.getcwd(), './4-quiver/', m_ctg_id)
        #ref_fasta = makePypeLocalFile(os.path.join(wd, '{ctg_id}_ref.fa'.format(ctg_id = ctg_id)))
        #read_bam = makePypeLocalFile(os.path.join(os.getcwd(), './4-quiver/reads/' '{ctg_id}.sam'.format(ctg_id = ctg_id)))
        ctg_type = ctg_types[ctg_id]
        cns_fasta = makePypeLocalFile(os.path.join(wd, 'cns-{ctg_id}.fasta.gz'.format(ctg_id=ctg_id)))
        cns_fastq = makePypeLocalFile(os.path.join(wd, 'cns-{ctg_id}.fastq.gz'.format(ctg_id=ctg_id)))
        job_done = makePypeLocalFile(os.path.join(wd, '{ctg_id}_{ctg_type}_quiver_done'.format(
            **locals())))
        if ctg_type == 'p':
            p_ctg_out.append((fn(cns_fasta), fn(cns_fastq)))
        elif ctg_type == 'h':
            h_ctg_out.append((fn(cns_fasta), fn(cns_fastq)))
        else:
            msg = 'Type is {!r}, not "p" or "h". Why are we running Quiver?'.format(ctg_type)
            raise Exception(msg)

        if os.path.exists(fn(read_bam)):  # TODO(CD): Ask Jason what we should do if missing SAM.
            parameters = {
                'job_uid': 'q-' + ctg_id,
                'ctg_id': ctg_id,
                'smrt_bin': smrt_bin,
                'sge_option': sge_option,
            }
            make_quiver_task = PypeTask(
                    inputs={
                        'ref_fasta': ref_fasta, 'read_bam': read_bam,
                        'scattered_quiver': scattered_quiver_plf,
                    },
                    outputs={
                        'cns_fasta': cns_fasta, 'cns_fastq': cns_fastq, 'job_done': job_done,
                    },
                    parameters=parameters,
                    )
            quiver_task = make_quiver_task(task_run_quiver)
            yield quiver_task
            job_done_plfs['{}'.format(ctg_id)] = job_done

    make_task = PypeTask(
        inputs=job_done_plfs,
        outputs={
            'gathered_quiver': gathered_quiver_plf,
        },
        parameters={},
    )
    yield make_task(task_gather_quiver)


def get_cns_zcat_task(
        gathered_quiver_plf,
        zcat_done_plf,
):
    cns_p_ctg_fasta_plf = makePypeLocalFile('4-quiver/cns_output/cns_p_ctg.fasta')
    cns_p_ctg_fastq_plf = makePypeLocalFile('4-quiver/cns_output/cns_p_ctg.fastq')
    cns_h_ctg_fasta_plf = makePypeLocalFile('4-quiver/cns_output/cns_h_ctg.fasta')
    cns_h_ctg_fastq_plf = makePypeLocalFile('4-quiver/cns_output/cns_h_ctg.fastq')
    make_task = PypeTask(
        inputs={
            'gathered_quiver': gathered_quiver_plf,
        },
        outputs={
            'cns_p_ctg_fasta': cns_p_ctg_fasta_plf,
            'cns_p_ctg_fastq': cns_p_ctg_fastq_plf,
            'cns_h_ctg_fasta': cns_h_ctg_fasta_plf,
            'cns_h_ctg_fastq': cns_h_ctg_fastq_plf,
            'job_done': zcat_done_plf,
        },
    )
    return make_task(task_cns_zcat)


def run_workflow(wf, config):
    abscwd = os.path.abspath('.')
    parameters = {
        'sge_option': config['sge_track_reads'],  # applies to select_reads task also, for now
        'max_n_open_files': config['max_n_open_files'],
        'topdir': os.getcwd(),
    }
    input_bam_fofn = config['input_bam_fofn']
    track_reads_h_done = './4-quiver/track_reads/track_reads_h_done'
    track_reads_rr2c = './4-quiver/track_reads/rawread_to_contigs'
    wf.addTask(gen_task(
        script=TASK_TRACK_READS_H_SCRIPT,
        inputs={
            'input_bam_fofn': input_bam_fofn,
            'hasm_done': './3-unzip/1-hasm/hasm_done',
        },
        outputs={
            'job_done': track_reads_h_done,
            'rawread_to_contigs': track_reads_rr2c,
        },
        parameters=parameters,
    ))

    read2ctg = './4-quiver/select_reads/read2ctg.msgpack'
    wf.addTask(gen_task(
        script=TASK_SELECT_READS_H_SCRIPT,
        inputs={
            # Some implicit inputs, plus these deps:
            'track_reads_h_done': track_reads_h_done,
            'input_bam_fofn': input_bam_fofn,
            #'rawread_to_contigs': track_reads_rr2c, # TODO: Check, and make explicit, maybe.
        },
        outputs={
            'read2ctg': read2ctg,
        },
        parameters=parameters,
    ))

    read2ctg_plf = makePypeLocalFile(read2ctg)
    input_bam_fofn_plf = makePypeLocalFile(input_bam_fofn)

    merged_fofn_plf = makePypeLocalFile('./4-quiver/merge_reads/merged.fofn')
    task = PypeTask(inputs={
        'input_bam_fofn': input_bam_fofn_plf,
        'read2ctg': read2ctg_plf},
        outputs={
        'merged_fofn': merged_fofn_plf},
        parameters=parameters,
    )(task_merge_reads)
    wf.addTask(task)

    scattered_segregate_plf = makePypeLocalFile('./4-quiver/segregate_scatter/scattered.json')
    task = PypeTask(
        inputs={
            'merged_fofn': merged_fofn_plf,
        },
        outputs={
            'scattered_segregate_json': scattered_segregate_plf,
        },
        parameters=parameters,
    )(task_segregate_scatter)
    wf.addTask(task)
    wf.refreshTargets()

    ctg2segregated_bamfn_plf = makePypeLocalFile('./4-quiver/segregate_gather/ctg2segregated_bamfn.msgpack')
    wf.addTasks(list(yield_segregate_bam_tasks(
        parameters, scattered_segregate_plf, ctg2segregated_bamfn_plf)))

    scattered_quiver_plf = makePypeLocalFile('4-quiver/quiver_scatter/scattered.json')
    parameters = {
        'config': config,
    }
    wf.addTask(get_scatter_quiver_task(
        parameters, ctg2segregated_bamfn_plf,
        scattered_quiver_plf,
        ))
    wf.refreshTargets()

    gathered_quiver_plf = makePypeLocalFile('4-quiver/cns_gather/gathered_quiver.json')

    wf.addTasks(list(yield_quiver_tasks(
        scattered_quiver_plf,
        gathered_quiver_plf)))

    zcat_done_plf = makePypeLocalFile('4-quiver/cns_output/job_done')

    wf.addTask(get_cns_zcat_task(
        gathered_quiver_plf,
        zcat_done_plf))

    wf.refreshTargets()
