from pypeflow.simple_pwatcher_bridge import (
    PypeLocalFile, makePypeLocalFile, fn,
    PypeTask,
)
import os


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


def get_phasing_tasks(phased_reads_file, bam_file, fasta_file, ctg_id, base_dir):
    vmap_file = makePypeLocalFile(os.path.join(base_dir, ctg_id, 'het_call', "variant_map"))
    vpos_file = makePypeLocalFile(os.path.join(base_dir, ctg_id, 'het_call', "variant_pos"))
    q_id_map_file = makePypeLocalFile(os.path.join(base_dir, ctg_id, 'het_call', "q_id_map"))
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


def task_phasing_readmap(self):
    job_done = fn(self.job_done)
    phased_reads_fn = fn(self.phased_reads)
    rid_to_phase_out_fn = fn(self.rid_to_phase_out)
    wd = self.parameters['wd']
    ctg_id = self.parameters['ctg_id']
    script_fn = os.path.join(wd, 'phasing_readmap_%s.sh' % (ctg_id))

    script = """\
set -vex
trap 'touch {job_done}.exit' EXIT
hostname
date
cd {wd}
python -m falcon_unzip.mains.phasing_readmap --the-ctg-id {ctg_id} --read-map-dir ../../../2-asm-falcon/read_maps --phased-reads {phased_reads_fn} >| {rid_to_phase_out_fn}.tmp
mv {rid_to_phase_out_fn}.tmp {rid_to_phase_out_fn}
date
touch {job_done}
""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def create_phasing_tasks(config, ctg_ids, all_ctg_out):
    """Report outputs via all_ctg_out.
    """
    for ctg_id in ctg_ids:
        # work-dir
        wd = os.path.join(os.getcwd(), './3-unzip/0-phasing/{ctg_id}/'.format(ctg_id=ctg_id))

        # inputs (by convention, for now)
        ref_fasta = makePypeLocalFile('./3-unzip/reads/{ctg_id}_ref.fa'.format(ctg_id=ctg_id))
        read_fasta = makePypeLocalFile('./3-unzip/reads/{ctg_id}_reads.fa'.format(ctg_id=ctg_id))

        blasr_dir = os.path.join(wd, 'blasr')
        ctg_aln_out = makePypeLocalFile(os.path.join(
            blasr_dir, '{ctg_id}_sorted.bam'.format(ctg_id=ctg_id)))

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
        for task in get_phasing_tasks(**kwds):
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
