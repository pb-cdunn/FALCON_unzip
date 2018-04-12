from __future__ import absolute_import
#from falcon_kit.pype import (wrap_gen_task as gen_task, gen_parallel_tasks, Dist)
from falcon_kit.pype import Dist
#from falcon_kit import pype_tasks
#from .pype import gen_task, gen_parallel_tasks
from falcon_kit.pype import (wrap_gen_task as gen_task, gen_parallel_tasks, Dist)
from .. import io
import logging
import os
import re

LOG = logging.getLogger(__name__)


# For now, in/outputs are in various directories, by convention, including '0-rawreads/m_*/*.msgpack'
TASK_TRACK_READS_H_SCRIPT = """\
python -m falcon_unzip.mains.get_read_hctg_map --base-dir={params.topdir} --output=read_to_contig_map
# formerly generated ./4-quiver/read_maps/read_to_contig_map

fc_rr_hctg_track.py --base-dir={params.topdir} --stream --read-to-contig-map=read_to_contig_map
# That writes into 0-rawreads/m_*/
# n_core is actually limited by number of files, but in theory we could use whole machine,
# Note: We also use a proc for LA4Falcon, so this is half.

abs_rawread_to_contigs=$(readlink -f {output.rawread_to_contigs}) #TODO: No readlink
cwd=$(pwd)
cd {params.topdir}
fc_rr_hctg_track2.exe --read-to-contig-map=${{cwd}}/read_to_contig_map --output=${{abs_rawread_to_contigs}}
cd -
ls -l {output.rawread_to_contigs}
"""


# For now, in/outputs are in various directories, by convention.
TASK_SELECT_READS_H_SCRIPT = """\
python -m falcon_unzip.mains.get_read2ctg --rawread-to-contigs={input.rawread_to_contigs} --output={output.read2ctg} {input.input_bam_fofn}
# I think this should be memory-constrained. Only 1 proc.
"""

# For now, in/outputs are in various directories, by convention.
TASK_MERGE_READS_SCRIPT = """\
rm -f {output.merged_fofn}
abs_output_merged_fofn=$(readlink -f {output.merged_fofn})
abs_input_input_bam_fofn=$(readlink -f {input.input_bam_fofn})
abs_input_read2ctg=$(readlink -f {input.read2ctg})
cd {params.topdir}
pwd
#fc_select_reads_from_bam.py --max-n-open-files={params.max_n_open_files} ${{abs_input_input_bam_fofn}}
python -m falcon_unzip.mains.bam_partition_and_merge --max-n-open-files={params.max_n_open_files} --read2ctg-fn=${{abs_input_read2ctg}} --merged-fn=${{abs_output_merged_fofn}} ${{abs_input_input_bam_fofn}}
# I think this should be memory-constrained. Only 1 proc, I think.
cd -
ls -l {output.merged_fofn}
"""

TASK_MAP_SEGREGATED_BAM_SCRIPT = """
python -m falcon_unzip.mains.bam_segregate_gather --gathered-fn={input.gathered} --ctg2segregated-bamfn-fn={output.ctg2segregated_bamfn}
"""

TASK_QUIVER_RUN_SCRIPT = """\
set -vex
trap 'touch {output.job_done}.exit' EXIT
hostname
date

samtools faidx {input.ref_fasta}
pbalign --tmpDir=/scratch/ --nproc=24 --minAccuracy=0.75 --minLength=50\
          --minAnchorSize=12 --maxDivergence=30 --concordant --algorithm=blasr\
          --algorithmOptions=--useQuality --maxHits=1 --hitPolicy=random --seed=1\
            {input.read_bam} {input.ref_fasta} aln-{params.ctg_id}.bam
#python -c 'import ConsensusCore2 as cc2; print cc2' # So quiver likely works.
(variantCaller --algorithm=arrow -x 5 -X 120 -q 20 -j 24 -r {input.ref_fasta} aln-{params.ctg_id}.bam\
            -o {output.cns_fasta} -o {output.cns_fastq}) || echo WARNING quiver failed. Maybe no reads for this block.
touch {output.cns_fasta}
touch {output.cns_fastq}
cp -f {input.ctg_type} {output.ctg_type_again}
date
touch {output.job_done}
"""


TASK_CNS_ZCAT_SCRIPT = """\
python -m falcon_unzip.mains.cns_zcat \
    --gathered-quiver-fn={input.gathered_quiver} \
    --cns-p-ctg-fasta-fn={output.cns_p_ctg_fasta} \
    --cns-p-ctg-fastq-fn={output.cns_p_ctg_fastq} \
    --cns-h-ctg-fasta-fn={output.cns_h_ctg_fasta} \
    --cns-h-ctg-fastq-fn={output.cns_h_ctg_fastq} \

touch {output.job_done}
"""

TASK_QUIVER_SPLIT_SCRIPT = """\
python -m falcon_unzip.mains.quiver_split --p-ctg-fasta-fn={input.p_ctg_fa} --h-ctg-fasta-fn={input.h_ctg_fa} --ctg2bamfn-fn={input.ctg2bamfn} --split-fn={output.split} --bash-template-fn={output.bash_template}
"""

TASK_SEPARATE_GATHERED_QUIVER_SCRIPT = """\
python -m falcon_unzip.mains.quiver_separate_gathered --gathered-fn={input.gathered} --output-fn={output.separated}
"""

TASK_SEGREGATE_SPLIT_SCRIPT = """
python -m falcon_unzip.mains.bam_segregate_split --merged-fofn-fn={input.merged_fofn} --split-fn={output.split} --bash-template-fn={output.bash_template}
"""

TASK_SEGREGATE_RUN_SCRIPT = """
python -m falcon_unzip.mains.bam_segregate --merged-bam-fn={input.merged_bam_fn} --segregated-bam-fns-fn={output.segregated_bam_fns}
"""
# ctg is encoded into each filepath within the output FOFN.
#   e.g. './4-quiver/segregate_scatter/segr001/000000F/000000F.bam'
# max_n_open_files = 300 # Ignored for now. Should not matter here.


def run_workflow(wf, config, rule_writer):
    default_njobs = int(config['job.defaults']['njobs'])
    #import pdb; pdb.set_trace()
    input_bam_fofn = os.path.relpath(config['Unzip']['input_bam_fofn']) # All input paths should be relative, for snakemake.
    track_reads_rr2c = './4-quiver/track-reads/rawread_to_contigs'
    wf.addTask(gen_task(
        script=TASK_TRACK_READS_H_SCRIPT,
        inputs={
            #'input_bam_fofn': input_bam_fofn,
            'hasm_done': './3-unzip/hasm_done',
        },
        outputs={
            'rawread_to_contigs': track_reads_rr2c,
        },
        parameters={},
        rule_writer=rule_writer,
        dist=Dist(NPROC=12, # guesstimate
            job_dict=config['job.step.unzip.track_reads'],
            use_tmpdir=False,
        ),
    ))

    read2ctg = './4-quiver/select-reads/read2ctg.msgpack'
    wf.addTask(gen_task(
        script=TASK_SELECT_READS_H_SCRIPT,
        inputs={
            # Some implicit inputs, plus these deps:
            'input_bam_fofn': input_bam_fofn,
            'rawread_to_contigs': track_reads_rr2c,
        },
        outputs={
            'read2ctg': read2ctg,
        },
        parameters={},
        rule_writer=rule_writer,
        dist=Dist(NPROC=4, MB=16, # actually NPROC=1, but our qsub jobs rarely report mem needs
            job_dict=config['job.step.unzip.track_reads'],
            use_tmpdir=False,
        ),
    ))

    #read2ctg_plf = makePypeLocalFile(read2ctg)
    #input_bam_fofn_plf = makePypeLocalFile(input_bam_fofn)

    merged_fofn = './4-quiver/merge-reads/merged.fofn'
    wf.addTask(gen_task(
        script=TASK_MERGE_READS_SCRIPT,
        inputs={
            'input_bam_fofn': input_bam_fofn,
            'read2ctg': read2ctg,
        },
        outputs={
            'merged_fofn': merged_fofn,
        },
        parameters={
            'max_n_open_files': config['max_n_open_files'],
        },
        rule_writer=rule_writer,
        dist=Dist(NPROC=4, MB=16, # actually NPROC=1, but our qsub jobs rarely report mem needs
            job_dict=config['job.step.unzip.track_reads'],
            use_tmpdir=False, # until we ensure the output uses non-tmpdir paths
        ),
    ))

    segr_all_units_fn ='./4-quiver/segregate-split/all-units-of-work.json'
    segr_run_bash_template_fn ='./4-quiver/segregate-split/bash-template.sh'
    wf.addTask(gen_task(
        script=TASK_SEGREGATE_SPLIT_SCRIPT,
        inputs=dict(
            merged_fofn=merged_fofn,
        ),
        outputs=dict(
            split=segr_all_units_fn,
            bash_template=segr_run_bash_template_fn,
        ),
        parameters={},
        rule_writer=rule_writer,
        dist=Dist(local=True),
    ))

    # Segregate reads from merged BAM files in parallel.
    # (If this were not done in Python, it could probably be in serial.)
    gathered_fn = './4-quiver/segregate-gathered/segregated-bam.json'
    gen_parallel_tasks(
        wf, rule_writer,
        segr_all_units_fn, gathered_fn,
        run_dict=dict(
            bash_template_fn=segr_run_bash_template_fn,
            script=TASK_SEGREGATE_RUN_SCRIPT,
            inputs={
                'units_of_work': './4-quiver/segregate-chunks/{segr}/some-units-of-work.json',
            },
            outputs={
                'results': './4-quiver/segregate-run/{segr}/segregated-bam-fn-list.json',
            },
            parameters={},
        ),
        dist=Dist(NPROC=4, MB=16, # actually NPROC=1, but our qsub jobs rarely report mem needs
            job_dict=config['job.step.unzip.track_reads'],
            use_tmpdir=False,
        ),
    )

    ctg2segregated_bamfn = './4-quiver/segregated-bam/ctg2segregated_bamfn.msgpack'
    # This is a separate task, consuming the output of the implicit gatherer.
    wf.addTask(gen_task(
        script=TASK_MAP_SEGREGATED_BAM_SCRIPT,
        inputs={
            'gathered': gathered_fn,
        },
        outputs={
            'ctg2segregated_bamfn': ctg2segregated_bamfn,
        },
        parameters={},
        rule_writer=rule_writer,
        dist=Dist(local=True),
    ))

    quiver_all_units_fn ='./4-quiver/quiver-split/all-units-of-work.json'
    quiver_run_bash_template_fn ='./4-quiver/quiver-split/bash-template.sh'
    wf.addTask(gen_task(
        script=TASK_QUIVER_SPLIT_SCRIPT,
        inputs={
                'p_ctg_fa': './3-unzip/all_p_ctg.fa',
                'h_ctg_fa': './3-unzip/all_h_ctg.fa',
                'ctg2bamfn': ctg2segregated_bamfn,
        },
        outputs={
                'split': quiver_all_units_fn,
                'bash_template': quiver_run_bash_template_fn,
        },
        parameters={},
        rule_writer=rule_writer,
        dist=Dist(local=True), # TODO: lots of fasta parsing, but we must not run in tmpdir
    ))

    unzip_quiver_njobs = int(config['job.step.unzip.quiver'].get('njobs', default_njobs))
    wf.max_jobs = unzip_quiver_njobs

    int_gathered_fn = '4-quiver/cns-gather/intermediate/int.gathered.json'
    gen_parallel_tasks(
        wf, rule_writer,
        quiver_all_units_fn, int_gathered_fn,
        run_dict=dict(
            bash_template_fn=quiver_run_bash_template_fn,
            script=TASK_QUIVER_RUN_SCRIPT,
            inputs={
                'units_of_work': './4-quiver/quiver-chunks/{ctg_id}/some-units-of-work.json',
            },
            outputs={
                'results': './4-quiver/quiver-run/{ctg_id}/results.json',
            },
            parameters={},
        ),
        dist=Dist(NPROC=24,
            job_dict=config['job.step.unzip.quiver'],
        ),
    )
    gathered_quiver = '4-quiver/cns-gather/gathered_quiver.json'
    wf.addTask(gen_task(
        script=TASK_SEPARATE_GATHERED_QUIVER_SCRIPT,
        inputs={
            'gathered': int_gathered_fn,
        },
        outputs={
            'separated': gathered_quiver,
        },
        parameters={},
        rule_writer=rule_writer,
        dist=Dist(local=True),
    ))

    wf.addTask(gen_task(
        script=TASK_CNS_ZCAT_SCRIPT,
        inputs={
            'gathered_quiver': gathered_quiver,
        },
        outputs={
            'cns_p_ctg_fasta': '4-quiver/cns-output/cns_p_ctg.fasta',
            'cns_p_ctg_fastq': '4-quiver/cns-output/cns_p_ctg.fastq',
            'cns_h_ctg_fasta': '4-quiver/cns-output/cns_h_ctg.fasta',
            'cns_h_ctg_fastq': '4-quiver/cns-output/cns_h_ctg.fastq',
            'job_done': '4-quiver/cns-output/job_done',
        },
        rule_writer=rule_writer,
        dist=Dist(NPROC=1),
    ))

    wf.refreshTargets()
