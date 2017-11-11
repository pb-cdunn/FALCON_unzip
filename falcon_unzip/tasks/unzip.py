from pypeflow.simple_pwatcher_bridge import (
    PypeLocalFile, makePypeLocalFile, fn,
    PypeTask,
)
from falcon_kit import pype_tasks
import logging
import os

LOG = logging.getLogger(__name__)


def create_tasks_read_to_contig_map(read_to_contig_map_plf):
    falcon_asm_done = makePypeLocalFile('./2-asm-falcon/falcon_asm_done')

    rawread_db = makePypeLocalFile('0-rawreads/raw_reads.db')
    rawread_ids = makePypeLocalFile('3-unzip/reads/dump_rawread_ids/rawread_ids')

    task = PypeTask(
        inputs={'rawread_db': rawread_db,
                'falcon_asm_done': falcon_asm_done,
        },
        outputs={'rawread_id_file': rawread_ids,
        },
    )
    yield task(pype_tasks.task_dump_rawread_ids)

    pread_db = makePypeLocalFile('1-preads_ovl/preads.db')
    pread_ids = makePypeLocalFile('3-unzip/reads/dump_pread_ids/pread_ids')

    task = PypeTask(
        inputs={'pread_db': pread_db,
                'falcon_asm_done': falcon_asm_done,
        },
        outputs={'pread_id_file': pread_ids,
        },
    )
    yield task(pype_tasks.task_dump_pread_ids)

    sg_edges_list = makePypeLocalFile('2-asm-falcon/sg_edges_list')
    utg_data = makePypeLocalFile('2-asm-falcon/utg_data')
    ctg_paths = makePypeLocalFile('2-asm-falcon/ctg_paths')

    inputs = {'rawread_id_file': rawread_ids,
              'pread_id_file': pread_ids,
              'sg_edges_list': sg_edges_list,
              'utg_data': utg_data,
              'ctg_paths': ctg_paths}
    task = PypeTask(
        inputs=inputs,
        outputs={'read_to_contig_map': read_to_contig_map_plf},
    )
    yield task(pype_tasks.task_generate_read_to_ctg_map)


def task_run_blasr(self):
    ctg_aln_out = fn(self.ctg_aln_out)
    ref_fasta = fn(self.ref_fasta)
    read_fasta = fn(self.read_fasta)

    job_uid = self.parameters['job_uid']
    ctg_id = self.parameters['ctg_id']

    script_fn = 'aln_{ctg_id}.sh'.format(ctg_id=ctg_id)
    script = """\
set -vex
hostname
date
time blasr {read_fasta} {ref_fasta} --noSplitSubreads --clipping subread\
 --hitPolicy randombest --randomSeed 42 --bestn 1 --minPctIdentity 70.0\
 --minMatch 12  --nproc 24 --bam --out tmp_aln.bam
#samtools view -bS tmp_aln.sam | samtools sort - {ctg_id}_sorted
samtools sort tmp_aln.bam -o {ctg_aln_out}
samtools index {ctg_aln_out}
rm tmp_aln.bam
date
""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def get_blasr_task(config, ctg_id, ctg_aln_out):
        ref_fasta = makePypeLocalFile('./3-unzip/reads/{ctg_id}_ref.fa'.format(ctg_id=ctg_id))
        read_fasta = makePypeLocalFile('./3-unzip/reads/{ctg_id}_reads.fa'.format(ctg_id=ctg_id))

        parameters = {'job_uid': 'aln-' + ctg_id, 'ctg_id': ctg_id,
                      'sge_option': config['sge_blasr_aln'],
                      }
        make_blasr_task = PypeTask(
                inputs={
                    'ref_fasta': ref_fasta,
                    'read_fasta': read_fasta,
                },
                outputs={
                    'ctg_aln_out': ctg_aln_out,
                },
                parameters=parameters,
                )
        blasr_task = make_blasr_task(task_run_blasr)
        return blasr_task


def task_make_het_call(self):
    bam_fn = fn(self.bam_file)
    fasta_fn = fn(self.fasta)
    vmap_fn = fn(self.vmap_file)
    vpos_fn = fn(self.vpos_file)
    q_id_map_fn = fn(self.q_id_map_file)
    ctg_id = self.parameters['ctg_id']
    script_fn = 'phasing_make_het_call.sh'
    script = """\
set -vex
hostname
pwd
date
python -m falcon_unzip.mains.phasing_make_het_call --bam {bam_fn} --fasta {fasta_fn} --ctg-id {ctg_id} --vmap={vmap_fn} --vpos={vpos_fn} --q-id-map={q_id_map_fn}
date
""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def task_generate_association_table(self):
    vmap_fn = fn(self.vmap_file)
    atable_fn = fn(self.atable_file)
    ctg_id = self.parameters['ctg_id']
    script_fn = 'phasing_generate_association_table.sh'
    script = """\
set -vex
hostname
pwd
date
python -m falcon_unzip.mains.phasing_generate_association_table --ctg-id {ctg_id} --vmap={vmap_fn} --atable={atable_fn}
date
""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def task_get_phased_blocks(self):
    vmap_fn = fn(self.vmap_file)
    atable_fn = fn(self.atable_file)
    p_variant_fn = fn(self.phased_variant_file)
    script_fn = 'phasing_get_phased_blocks.sh'
    script = """\
set -vex
hostname
pwd
date
python -m falcon_unzip.mains.phasing_get_phased_blocks --vmap={vmap_fn} --atable={atable_fn} --p-variant={p_variant_fn}
date
""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def task_get_phased_reads(self):
    phased_reads_fn = fn(self.phased_reads_file)
    q_id_map_fn = fn(self.q_id_map_file)
    vmap_fn = fn(self.vmap_file)
    p_variant_fn = fn(self.phased_variant_file)
    ctg_id = self.parameters['ctg_id']
    script_fn = 'phasing_get_phased_reads.sh'
    script = """\
set -vex
hostname
pwd
date
python -m falcon_unzip.mains.phasing_get_phased_reads --ctg-id={ctg_id} --vmap={vmap_fn} --p-variant={p_variant_fn} --q-id-map={q_id_map_fn} --phased-reads={phased_reads_fn}
date
""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def task_phasing(self):
    ref_fasta = fn(self.ref_fasta)
    aln_bam = fn(self.aln_bam)

    job_done = fn(self.job_done)

    job_uid = self.parameters['job_uid']
    wd = self.parameters['wd']
    ctg_id = self.parameters['ctg_id']

    config = self.parameters['config']
    smrt_bin = config['smrt_bin']

    script_fn = os.path.join(wd, 'p_%s.sh' % (ctg_id))

    script = """\
set -vex
trap 'touch {job_done}.exit' EXIT
hostname
date
cd {wd}
mkdir -p phasing_subworkflow
cd phasing_subworkflow
python -m falcon_unzip.mains.phasing --bam {aln_bam} --fasta {ref_fasta} --ctg_id {ctg_id} --base_dir ../..
date
cd ..
touch {job_done}
""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def yield_phasing_tasks(phased_reads_file, bam_file, fasta_file, ctg_id, base_dir):
    vmap_file = makePypeLocalFile(os.path.join(base_dir, ctg_id, 'het_call', "variant_map"))
    vpos_file = makePypeLocalFile(os.path.join(base_dir, ctg_id, 'het_call', "variant_pos"))
    q_id_map_file = makePypeLocalFile(os.path.join(base_dir, ctg_id, 'het_call', "q_id_map.msgpack"))
    parameters = {}
    parameters["ctg_id"] = ctg_id
    parameters["base_dir"] = base_dir

    make_het_call_task = PypeTask(
        inputs={
            "bam_file": bam_file,
            "fasta": fasta_file,
        },
        outputs={"vmap_file": vmap_file, "vpos_file": vpos_file, "q_id_map_file": q_id_map_file},
        parameters=parameters,
    )(task_make_het_call)

    yield make_het_call_task

    atable_file = makePypeLocalFile(os.path.join(base_dir, ctg_id, 'g_atable', "atable"))
    parameters = {}
    parameters["ctg_id"] = ctg_id
    parameters["base_dir"] = base_dir
    generate_association_table_task = PypeTask(inputs={"vmap_file": vmap_file},
                                               outputs={"atable_file": atable_file},
                                               parameters=parameters,
                                               )(task_generate_association_table)

    yield generate_association_table_task

    phased_variant_file = makePypeLocalFile(os.path.join(base_dir, ctg_id, 'get_phased_blocks', "phased_variants"))
    get_phased_blocks_task = PypeTask(inputs={"vmap_file": vmap_file, "atable_file": atable_file},
                                      outputs={"phased_variant_file": phased_variant_file},
                                      )(task_get_phased_blocks)
    yield get_phased_blocks_task

    get_phased_reads_task = PypeTask(inputs={"vmap_file": vmap_file,
                                             "q_id_map_file": q_id_map_file,
                                             "phased_variant_file": phased_variant_file},
                                     outputs={"phased_reads_file": phased_reads_file},
                                     parameters={"ctg_id": ctg_id},
                                     )(task_get_phased_reads)
    yield get_phased_reads_task


def task_get_rid_to_phase_all(self):
    rid_to_phase_all = fn(self.rid_to_phase_all)
    inputs_fn = [fn(f) for f in self.inputs.values()]
    inputs_fn.sort()
    input = ' '.join(i for i in inputs_fn)
    LOG.info('Generate {!r} from {!r}'.format(
        rid_to_phase_all, input))
    script = """
rm -f {rid_to_phase_all}
for fn in {input}; do
  cat $fn >> {rid_to_phase_all}
done
""".format(**locals())
    script_fn = 'gather_rid_to_phase.sh'
    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def task_phasing_readmap(self):
    # TODO: read-map-dir/* as inputs
    job_done = fn(self.job_done)
    phased_reads_fn = fn(self.phased_reads)
    rid_to_phase_out_fn = fn(self.rid_to_phase_out)

    ctg_id = self.parameters['ctg_id']

    script_fn = 'phasing_readmap_{}.sh'.format(ctg_id)
    script = """\
set -vex
trap 'touch {job_done}.exit' EXIT
hostname
date
python -m falcon_unzip.mains.phasing_readmap --the-ctg-id {ctg_id} --read-map-dir ../../reads --phased-reads {phased_reads_fn} >| {rid_to_phase_out_fn}.tmp
mv {rid_to_phase_out_fn}.tmp {rid_to_phase_out_fn}
date
touch {job_done}
""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def create_phasing_tasks(config, ctg_list_file, rid_to_phase_all):
    """
    Gathered input is ctg_list_file.
    Gathered output is rid_to_phase_all.
    """
    ctg_ids = []
    with open(fn(ctg_list_file)) as f:
        for row in f:
            row = row.strip()
            ctg_ids.append(row)

    all_ctg_out = dict()

    for ctg_id in ctg_ids:
        # work-dir
        wd = './3-unzip/0-phasing/{ctg_id}/'.format(ctg_id=ctg_id)
        ctg_aln_out = makePypeLocalFile(
            '{wd}/blasr/{ctg_id}_sorted.bam'.format(wd=wd, ctg_id=ctg_id))

        yield get_blasr_task(config, ctg_id, ctg_aln_out)

        # inputs of basic phasing tasks
        ref_fasta = makePypeLocalFile('./3-unzip/reads/{ctg_id}_ref.fa'.format(ctg_id=ctg_id))
        read_fasta = makePypeLocalFile('./3-unzip/reads/{ctg_id}_reads.fa'.format(ctg_id=ctg_id))

        blasr_dir = os.path.join(wd, 'blasr')

        # output of basic phasing tasks
        phased_reads_file = makePypeLocalFile(os.path.join(
            wd, 'get_phased_reads', 'phased_reads'))

        kwds = dict(
            phased_reads_file=phased_reads_file,
            bam_file=ctg_aln_out,
            fasta_file=ref_fasta,
            ctg_id=ctg_id,
            base_dir='./3-unzip/0-phasing'
        )
        for task in yield_phasing_tasks(**kwds):
            yield task

        # final outputs
        job_done = makePypeLocalFile(os.path.join(
            wd, 'phasing_readmap_{ctg_id}_done'.format(ctg_id=ctg_id)))
        rid_to_phase_out = makePypeLocalFile(os.path.join(
            wd, 'rid_to_phase.{ctg_id}'.format(ctg_id=ctg_id)))
        all_ctg_out['r2p.{ctg_id}'.format(ctg_id=ctg_id)] = rid_to_phase_out

        parameters = {
            'wd': wd,
            'ctg_id': ctg_id,
        }
        make_task = PypeTask(
            inputs={'phased_reads': phased_reads_file,
                    },
            outputs={'job_done': job_done,
                     'rid_to_phase_out': rid_to_phase_out,
                     },
            parameters=parameters,
        )
        task = make_task(task_phasing_readmap)
        yield task

    task = PypeTask(
            inputs=all_ctg_out,
            outputs={
                'rid_to_phase_all': rid_to_phase_all,
            },
            )(task_get_rid_to_phase_all)
    yield task
