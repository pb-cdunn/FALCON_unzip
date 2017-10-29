from pypeflow.simple_pwatcher_bridge import (
    PypeLocalFile, makePypeLocalFile, fn,
    PypeTask,
)
import os


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
cd {wd}
hostname
date
cd {wd}
python -m falcon_unzip.phasing --bam {aln_bam} --fasta {ref_fasta} --ctg_id {ctg_id} --base_dir .. --samtools {samtools}
python -m falcon_unzip.phasing_readmap --ctg_id {ctg_id} --read_map_dir ../../../2-asm-falcon/read_maps --phased_reads phased_reads
date
touch {job_done}
""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


def create_phasing_tasks(config, ctg_ids, all_ctg_out):
    for ctg_id in ctg_ids:
        # inputs
        ref_fasta = makePypeLocalFile('./3-unzip/reads/{ctg_id}_ref.fa'.format(ctg_id=ctg_id))
        read_fasta = makePypeLocalFile('./3-unzip/reads/{ctg_id}_reads.fa'.format(ctg_id=ctg_id))

        # outputs
        wd = os.path.join(os.getcwd(), './3-unzip/0-phasing/{ctg_id}/'.format(ctg_id=ctg_id))

        blasr_dir = os.path.join(wd, 'blasr')
        ctg_aln_out = makePypeLocalFile(os.path.join(blasr_dir, '{ctg_id}_sorted.bam'.format(ctg_id=ctg_id)))

        phasing_dir = os.path.join(wd, 'phasing')
        job_done = makePypeLocalFile(os.path.join(phasing_dir, 'p_{ctg_id}_done'.format(ctg_id=ctg_id)))
        rid_to_phase_out = makePypeLocalFile(os.path.join(
            wd, 'rid_to_phase.{ctg_id}'.format(ctg_id=ctg_id)))  # TODO: ???
        all_ctg_out['r2p.{ctg_id}'.format(ctg_id=ctg_id)] = rid_to_phase_out  # implicit output?

        parameters = {'job_uid': 'ha-' + ctg_id, 'wd': wd, 'config': config, 'ctg_id': ctg_id,
                      'sge_option': config['sge_phasing'],
                      }
        make_phasing_task = PypeTask(inputs={'ref_fasta': ref_fasta, 'aln_bam': ctg_aln_out},
                                     outputs={'job_done': job_done},
                                     parameters=parameters,
                                     )
        phasing_task = make_phasing_task(task_phasing)
        return [phasing_task]
