from pypeflow.simple_pwatcher_bridge import (
    makePypeLocalFile, fn,
    PypeTask,
)
from falcon_kit.FastaReader import FastaReader
from .. import io
import json
import logging
import os

LOG = logging.getLogger(__name__)


def task_track_reads_h(self):
    input_bam_fofn = fn(self.input_bam_fofn)
    job_done = fn(self.job_done)
    topdir = os.path.relpath(self.parameters['topdir'])
    basedir = os.path.reldir(topdir)
    reldir = os.path.relpath('.', topdir)
    script_fn = 'track_reads_h.sh'

    # For now, in/outputs are in various directories, by convention, including '0-rawreads/m_*/*.msgpack'
    script = """\
set -vex
trap 'touch {job_done}.exit' EXIT
hostname
date

python -m falcon_unzip.mains.get_read_hctg_map --base-dir={basedir} --output=read_to_contig_map
# formerly generated ./4-quiver/read_maps/read_to_contig_map

rm -f ./3-unzip/reads/dump_rawread_ids/rawread_to_contigs

fc_rr_hctg_track.py --base-dir={basedir} --stream

cd {topdir}
fc_rr_hctg_track2.exe --output={reldir}/rawread_to_contigs
cd -

date
ls -l rawread_to_contigs
ls -lhH rawread_to_contigs
touch {job_done}
""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def task_select_reads_h(self):
    read2ctg_fn = fn(self.read2ctg)
    input_bam_fofn = fn(self.input_bam_fofn)
    topdir = os.path.relpath(self.parameters['topdir'])
    script_fn = 'select_reads_h.sh'

    # For now, in/outputs are in various directories, by convention.
    script = """\
set -vex
hostname
date

cd {topdir}
pwd
python -m falcon_unzip.mains.get_read2ctg --output={read2ctg_fn} {input_bam_fofn}

date
cd -
""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def task_merge_reads(self):
    merged_fofn_fn = fn(self.merged_fofn)
    read2ctg_fn = fn(self.read2ctg)
    input_bam_fofn = fn(self.input_bam_fofn)
    max_n_open_files = self.parameters['max_n_open_files']
    topdir = os.path.relpath(self.parameters['topdir'])
    script_fn = 'merge_reads.sh'

    # For now, in/outputs are in various directories, by convention.
    script = """\
set -vex
#trap 'touch {merged_fofn_fn}.exit' EXIT
hostname
date

cd {topdir}
#fc_select_reads_from_bam.py --max-n-open-files={max_n_open_files} {input_bam_fofn}
pwd
python -m falcon_unzip.mains.bam_partition_and_merge --max-n-open-files={max_n_open_files} --read2ctg-fn={read2ctg_fn} --merged-fn={merged_fofn_fn} {input_bam_fofn}

date
cd -
# Expect {merged_fofn_fn}
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
    gathered_p_ctg_fn = fn(self.gathered_p_ctg)
    gathered_h_ctg_fn = fn(self.gathered_h_ctg)

    cns_p_ctg_fasta_fn = fn(self.cns_p_ctg_fasta)
    cns_p_ctg_fastq_fn = fn(self.cns_p_ctg_fastq)
    cns_h_ctg_fasta_fn = fn(self.cns_h_ctg_fasta)
    cns_h_ctg_fastq_fn = fn(self.cns_h_ctg_fastq)
    job_done_fn = fn(self.job_done)

    io.rm(cns_p_ctg_fasta_fn)
    io.touch(cns_p_ctg_fasta_fn)
    io.rm(cns_p_ctg_fastq_fn)
    io.touch(cns_p_ctg_fastq_fn)
    with open(gathered_p_ctg_fn) as ifs:
        for line in ifs:
            cns_fasta_fn, cns_fastq_fn = line.split()
            io.syscall('zcat {cns_fasta_fn} >> {cns_p_ctg_fasta_fn}'.format(**locals()))
            io.syscall('zcat {cns_fastq_fn} >> {cns_p_ctg_fastq_fn}'.format(**locals()))

    # comment out this for now for recovering purpose
    # with open(gathered_p_ctg_fn) as ifs:
    #    for line in ifs:
    #        cns_fasta_fn, cns_fastq_fn = line.split()
    #        io.rm(cns_fasta_fn)
    #        io.rm(cns_fasta_fn)

    io.rm(cns_h_ctg_fasta_fn)
    io.touch(cns_h_ctg_fasta_fn)
    io.rm(cns_h_ctg_fastq_fn)
    io.touch(cns_h_ctg_fastq_fn)
    with open(gathered_h_ctg_fn) as ifs:
        for line in ifs:
            cns_fasta_fn, cns_fastq_fn = line.split()
            io.syscall('zcat {cns_fasta_fn} >> {cns_h_ctg_fasta_fn}'.format(**locals()))
            io.syscall('zcat {cns_fastq_fn} >> {cns_h_ctg_fastq_fn}'.format(**locals()))

    # comment out this for now for recovering purpose
    # with open(gathered_h_ctg_fn) as ifs:
    #    for line in ifs:
    #        cns_fasta_fn, cns_fastq_fn = line.split()
    #        io.rm(cns_fasta_fn)
    #        io.rm(cns_fasta_fn)

    io.touch(job_done_fn)


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
    job_done_fn = fn(self.job_done)
    io.touch(job_done_fn)


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
    merged_bamfn_fn = self.merged_bamfn
    segregated_bam_fofn_fn = self.segregated_bam_fofn

    script = """
python -m falcon_unzip.mains.bam_segregate --merged-fn={merged_bamfn_fn} --output-fn={segregated_bam_fofn_fn}
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


def create_quiver_jobs(wf, scattered_quiver_plf):
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
        cns_fasta = makePypeLocalFile(os.path.join(wd, 'cns-{ctg_id}.fasta.gz'.format(ctg_id=ctg_id)))
        cns_fastq = makePypeLocalFile(os.path.join(wd, 'cns-{ctg_id}.fastq.gz'.format(ctg_id=ctg_id)))
        job_done = makePypeLocalFile(os.path.join(wd, '{ctg_id}_quiver_done'.format(ctg_id=ctg_id)))

        if os.path.exists(fn(read_bam)):  # TODO(CD): Ask Jason what we should do if missing SAM.
            if ctg_types[ctg_id] == 'p':
                p_ctg_out.append((fn(cns_fasta), fn(cns_fastq)))
            elif ctg_types[ctg_id] == 'h':
                h_ctg_out.append((fn(cns_fasta), fn(cns_fastq)))
            else:
                LOG.warning('Type is {!r}, not "p" or "h". Why are we running Quiver?'.format(ctg_types[ctg_id]))
            parameters = {
                'job_uid': 'q-' + ctg_id,
                'ctg_id': ctg_id,
                'smrt_bin': smrt_bin,
                'sge_option': sge_option,
            }
            make_quiver_task = PypeTask(inputs={'ref_fasta': ref_fasta, 'read_bam': read_bam,
                                                'scattered_quiver': scattered_quiver_plf,
                                                },
                                        outputs={'cns_fasta': cns_fasta, 'cns_fastq': cns_fastq, 'job_done': job_done},
                                        parameters=parameters,
                                        )
            quiver_task = make_quiver_task(task_run_quiver)
            wf.addTask(quiver_task)
            job_done_plfs['{}'.format(ctg_id)] = job_done
    return p_ctg_out, h_ctg_out, job_done_plfs


def create_segregate_jobs(wf, parameters, scattered_segregate_plf):
    jn2segregated_bam_fofn = dict()  # job_name -> FOFN_plf

    #cwd = os.getcwd()
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
        wf.addTask(make_task(task_run_segregate))
        jn2segregated_bam_fofn[job_name] = segregated_bam_fofn_plf
    return jn2segregated_bam_fofn
