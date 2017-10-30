from pypeflow.simple_pwatcher_bridge import (
    PypeLocalFile, makePypeLocalFile, fn,
    PypeTask,
)
import os

from falcon_unzip import phasing

def task_make_het_call(self):
    bam_fn = fn(self.bam_file)
    fasta_fn = fn(self.fasta)
    vmap_fn = fn(self.vmap_file)
    vpos_fn = fn(self.vpos_file)
    q_id_map_fn = fn(self.q_id_map_file)
    kwds = dict(
            bam_fn=bam_fn,
            fasta_fn=fasta_fn,
            vmap_fn=vmap_fn,
            vpos_fn=vpos_fn,
            q_id_map_fn=q_id_map_fn,
            parameters=self.parameters)
    return phasing.make_het_call(**kwds)


def task_generate_association_table(self):
    vmap_fn = fn(self.vmap_file)
    atable_fn = fn(self.atable_file)
    kwds = dict(
            vmap_fn=vmap_fn,
            atable_fn=atable_fn,
            parameters=self.parameters)
    return phasing.generate_association_table(**kwds)


def task_get_phased_blocks(self):
    vmap_fn = fn(self.vmap_file)
    atable_fn = fn(self.atable_file)
    p_variant_fn = fn(self.phased_variant_file)
    kwds = dict(
            vmap_fn=vmap_fn,
            atable_fn=atable_fn,
            p_variant_fn=p_variant_fn)
    return phasing.get_phased_blocks(**kwds)


def task_get_phased_reads(self):
    phased_read_fn = fn(self.phased_read_file)
    q_id_map_fn = fn(self.q_id_map_file)
    vmap_fn = fn(self.vmap_file)
    p_variant_fn = fn(self.phased_variant_file)
    parameters = self.parameters
    kwds = dict(
            phased_read_fn=phased_read_fn,
            q_id_map_fn=q_id_map_fn,
            vmap_fn=vmap_fn,
            p_variant_fn=p_variant_fn,
            parameters=self.parameters)
    return phasing.get_phased_reads(**kwds)

def task_phasing(self):
    ref_fasta = fn(self.ref_fasta)
    aln_bam = fn(self.aln_bam)

    job_done = fn(self.job_done)

    job_uid = self.parameters['job_uid']
    wd = self.parameters['wd']
    ctg_id = self.parameters['ctg_id']

    config = self.parameters['config']
    smrt_bin = config['smrt_bin']
    samtools = os.path.join(smrt_bin, 'samtools')

    script_fn = os.path.join(wd, 'p_%s.sh' % (ctg_id))

    script = """\
set -vex
trap 'touch {job_done}.exit' EXIT
hostname
date
cd {wd}
mkdir -p phasing_subworkflow
cd phasing_subworkflow
python -m falcon_unzip.mains.phasing --bam {aln_bam} --fasta {ref_fasta} --ctg_id {ctg_id} --base_dir ../.. --samtools {samtools}
date
cd ..
touch {job_done}
""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def get_phasing_tasks(bam_fn, fasta_fn, ctg_id, base_dir, samtools):
    bam_file = makePypeLocalFile(bam_fn)
    fasta_file = makePypeLocalFile(fasta_fn)
    vmap_file = makePypeLocalFile(os.path.join(base_dir, ctg_id, 'het_call', "variant_map"))
    vpos_file = makePypeLocalFile(os.path.join(base_dir, ctg_id, 'het_call', "variant_pos"))
    q_id_map_file = makePypeLocalFile(os.path.join(base_dir, ctg_id, 'het_call', "q_id_map"))
    parameters = {}
    parameters["ctg_id"] = ctg_id
    parameters["base_dir"] = base_dir
    parameters["samtools"] = samtools

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

    phased_read_file = makePypeLocalFile(os.path.join(base_dir, ctg_id, "phased_reads"))
    get_phased_reads_task = PypeTask(inputs={"vmap_file": vmap_file,
                                             "q_id_map_file": q_id_map_file,
                                             "phased_variant_file": phased_variant_file},
                                     outputs={"phased_read_file": phased_read_file},
                                     parameters={"ctg_id": ctg_id},
                                     )(task_get_phased_reads)
    yield get_phased_reads_task


def task_phasing_readmap(self):
    job_done = fn(self.job_done)
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
python -m falcon_unzip.mains.phasing_readmap --the-ctg-id {ctg_id} --read-map-dir ../../../2-asm-falcon/read_maps --phased-reads phased_reads >| {rid_to_phase_out_fn}.tmp
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

        # outputs
        phasing_dir = os.path.join(wd, 'phasing')
        phasing_job_done = makePypeLocalFile(os.path.join(
            phasing_dir, 'p_{ctg_id}_done'.format(ctg_id=ctg_id)))
        rid_to_phase_out = makePypeLocalFile(os.path.join(
            wd, 'rid_to_phase.{ctg_id}'.format(ctg_id=ctg_id)))
        all_ctg_out['r2p.{ctg_id}'.format(ctg_id=ctg_id)] = rid_to_phase_out

        parameters = {'job_uid': 'ha-' + ctg_id, 'wd': wd, 'config': config, 'ctg_id': ctg_id,
                      'sge_option': config['sge_phasing'],
                      }
        make_task = PypeTask(
                inputs={'ref_fasta': ref_fasta, 'aln_bam': ctg_aln_out},
                outputs={'job_done': phasing_job_done},
                parameters=parameters,
        )
        task = make_task(task_phasing)
        yield task

        job_done = makePypeLocalFile(os.path.join(wd, 'phasing_readmap_{ctg_id}_done'.format(ctg_id=ctg_id)))
        make_task = PypeTask(
                inputs={'phasing_job_done': phasing_job_done,
                },
                outputs={'job_done': job_done,
                         'rid_to_phase_out': rid_to_phase_out,
                },
                parameters=parameters,
        )
        task = make_task(task_phasing_readmap)
        yield task
