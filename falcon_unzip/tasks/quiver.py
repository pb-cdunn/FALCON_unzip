from __future__ import absolute_import
from pypeflow.simple_pwatcher_bridge import (
    makePypeLocalFile, fn,
    PypeTask,
)
from .pype import (gen_task, gen_parallel_tasks)
from .. import io
import json
import logging
import os
import re

LOG = logging.getLogger(__name__)


# For now, in/outputs are in various directories, by convention, including '0-rawreads/m_*/*.msgpack'
TASK_TRACK_READS_H_SCRIPT = """\
python -m falcon_unzip.mains.get_read_hctg_map --base-dir={params.topdir} --output=read_to_contig_map
# formerly generated ./4-quiver/read_maps/read_to_contig_map

fc_rr_hctg_track.py --base-dir={params.topdir} --stream
# That writes into 0-rawreads/m_*/

abs_rawread_to_contigs=$(readlink -f {output.rawread_to_contigs}) #TODO: No readlink
cwd=$(pwd)
cd {params.topdir}
fc_rr_hctg_track2.exe --read-to-contig-map=${{cwd}}/read_to_contig_map --output=${{abs_rawread_to_contigs}}
cd -
ls -l {output.rawread_to_contigs}
"""


# For now, in/outputs are in various directories, by convention.
TASK_SELECT_READS_H_SCRIPT = """\
python -m falcon_unzip.mains.get_read2ctg --output={output.read2ctg} {input.input_bam_fofn}
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
cd -
ls -l {output.merged_fofn}
"""

TASK_MAP_SEGREGATED_BAM_SCRIPT = """
python -m falcon_unzip.mains.bam_segregate_gather --gathered-fn={input.gathered} --ctg2segregated-bamfn-fn={output.ctg2segregated_bamfn}
"""

TASK_RUN_QUIVER_SCRIPT = """\
set -vex
trap 'touch {output.job_done}.exit' EXIT
hostname
date

samtools faidx {input.ref_fasta}
pbalign --tmpDir=/localdisk/scratch/ --nproc=24 --minAccuracy=0.75 --minLength=50\
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

TASK_SCATTER_QUIVER_SCRIPT = """\
python -m falcon_unzip.mains.quiver_scatter --p-ctg-fasta-fn={input.p_ctg_fa} --h-ctg-fasta-fn={input.h_ctg_fa} --ctg2bamfn-fn={input.ctg2bamfn} --scattered-fn={output.scattered}
"""

TASK_SEPARATE_GATHERED_QUIVER_SCRIPT = """\
python -m falcon_unzip.mains.quiver_separate_gathered --gathered-fn={input.gathered} --output-fn={output.separated}
"""

TASK_SEGREGATE_SCATTER_SCRIPT = """
python -m falcon_unzip.mains.bam_segregate_scatter --merged-fofn-fn={input.merged_fofn} --scattered-fn={output.scattered}
"""

TASK_RUN_SEGREGATE_SCRIPT = """
python -m falcon_unzip.mains.bam_segregate --merged-fn={input.merged_bamfn} --output-fn={output.segregated_bam_fofn}
"""
# ctg is encoded into each filepath within the output FOFN.
#   e.g. './4-quiver/segregate_scatter/segr001/000000F/000000F.bam'
# max_n_open_files = 300 # Ignored for now. Should not matter here.


def run_workflow(wf, config, rule_writer):
    #import pdb; pdb.set_trace()
    parameters = {
        'sge_option': config['sge_track_reads'],  # applies to select_reads task also, for now
        'max_n_open_files': config['max_n_open_files'],
        'topdir': os.getcwd(),
    }
    input_bam_fofn = os.path.relpath(config['input_bam_fofn']) # All input paths should be relative, for snakemake.
    track_reads_rr2c = './4-quiver/track_reads/rawread_to_contigs'
    wf.addTask(gen_task(
        script=TASK_TRACK_READS_H_SCRIPT,
        inputs={
            #'input_bam_fofn': input_bam_fofn,
            'hasm_done': './3-unzip/1-hasm/hasm_done',
        },
        outputs={
            'rawread_to_contigs': track_reads_rr2c,
        },
        parameters=parameters,
        rule_writer=rule_writer,
    ))

    read2ctg = './4-quiver/select_reads/read2ctg.msgpack'
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
        parameters=parameters,
        rule_writer=rule_writer,
    ))

    #read2ctg_plf = makePypeLocalFile(read2ctg)
    #input_bam_fofn_plf = makePypeLocalFile(input_bam_fofn)

    merged_fofn = './4-quiver/merge_reads/merged.fofn'
    wf.addTask(gen_task(
        script=TASK_MERGE_READS_SCRIPT,
        inputs={
            'input_bam_fofn': input_bam_fofn,
            'read2ctg': read2ctg,
        },
        outputs={
            'merged_fofn': merged_fofn,
        },
        parameters=parameters,
        rule_writer=rule_writer,
    ))

    scattered ='./4-quiver/segregate_scatter/scattered.json'
    wf.addTask(gen_task(
        script=TASK_SEGREGATE_SCATTER_SCRIPT,
        inputs=dict(
            merged_fofn=merged_fofn,
        ),
        outputs=dict(
            scattered=scattered,
        ),
        parameters={},
        rule_writer=rule_writer,
    ))

    # Segregate reads from merged BAM files in parallel.
    # (If this were not done in Python, it could probably be in serial.)
    gathered = './4-quiver/segregate_gather/segregated_bam.json'
    gen_parallel_tasks(
        wf, rule_writer,
        scattered, gathered,
        run_dict=dict(
            script=TASK_RUN_SEGREGATE_SCRIPT,
            inputs={
                'merged_bamfn': './4-quiver/merge_reads/{segr}/merged.bam',
            },
            outputs={
                'segregated_bam_fofn': './4-quiver/segregate_scatter/{segr}/segregated_bam.fofn',
            },
            parameters=parameters,
        ),
    )

    ctg2segregated_bamfn = './4-quiver/segregate_bam/ctg2segregated_bamfn.msgpack'
    # This is a separate task, consuming the output of the implicit gatherer.
    wf.addTask(gen_task(
        script=TASK_MAP_SEGREGATED_BAM_SCRIPT,
        inputs={
            'gathered': gathered,
        },
        outputs={
            'ctg2segregated_bamfn': ctg2segregated_bamfn,
        },
        parameters=parameters,
        rule_writer=rule_writer,
    ))

    scattered_quiver = '4-quiver/quiver_scatter/scattered.json'
    wf.addTask(gen_task(
        script=TASK_SCATTER_QUIVER_SCRIPT,
        inputs={
                'p_ctg_fa': './3-unzip/all_p_ctg.fa',
                'h_ctg_fa': './3-unzip/all_h_ctg.fa',
                'ctg2bamfn': ctg2segregated_bamfn,
        },
        outputs={
                'scattered': scattered_quiver,
        },
        parameters={},
        rule_writer=rule_writer,
    ))

    int_gathered_fn = '4-quiver/cns_gather/intermediate/int.gathered.json'
    gen_parallel_tasks(
        wf, rule_writer,
        scattered_quiver, int_gathered_fn,
        run_dict=dict(
            script=TASK_RUN_QUIVER_SCRIPT,
            inputs={
                'read_bam': '4-quiver/segregate_scatter/segregated/{ctg_id}/reads.bam',
                'ref_fasta': '4-quiver/quiver_scatter/refs/{ctg_id}/ref.fa',
                'ctg_type': '4-quiver/quiver_scatter/refs/{ctg_id}/ctg_type',
            },
            outputs={
                'cns_fasta': '4-quiver/quiver_run/{ctg_id}/cns.fasta.gz',
                'cns_fastq': '4-quiver/quiver_run/{ctg_id}/cns.fastq.gz',
                'ctg_type_again': '4-quiver/quiver_run/{ctg_id}/ctg_type',
                'job_done': '4-quiver/quiver_run/{ctg_id}/quiver_done',
            },
            parameters=parameters, # expanded wildcards are added implicitly
            # TODO(CD): sge_quiver
        ),
    )
    gathered_quiver = '4-quiver/cns_gather/gathered_quiver.json'
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
    ))

    wf.addTask(gen_task(
        script=TASK_CNS_ZCAT_SCRIPT,
        inputs={
            'gathered_quiver': gathered_quiver,
        },
        outputs={
            'cns_p_ctg_fasta': '4-quiver/cns_output/cns_p_ctg.fasta',
            'cns_p_ctg_fastq': '4-quiver/cns_output/cns_p_ctg.fastq',
            'cns_h_ctg_fasta': '4-quiver/cns_output/cns_h_ctg.fasta',
            'cns_h_ctg_fastq': '4-quiver/cns_output/cns_h_ctg.fastq',
            'job_done': '4-quiver/cns_output/job_done',
        },
        rule_writer=rule_writer,
    ))

    wf.refreshTargets()
