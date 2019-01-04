from __future__ import absolute_import
from falcon_kit.pype import Dist
from falcon_kit import pype_tasks
from falcon_kit.pype import (wrap_gen_task as gen_task, gen_parallel_tasks, Dist)
from .. import io
import logging
import os

LOG = logging.getLogger(__name__)


TASK_TRACK_READS_SCRIPT = """\
nproc={params.pypeflow_nproc}
ncore=$(expr $(expr $(expr $nproc - 1 ) / 2) + 1)

# Also require read_to_contig_map.
python -m falcon_unzip.mains.rr_ctg_track --n-core=$ncore --db={input.r_db} --las-fofn={input.r_las_fofn} --output={output.rawread_to_contigs}
python -m falcon_unzip.mains.pr_ctg_track --n-core=$ncore --db={input.p_db} --las-fofn={input.p_las_fofn} --output={output.pread_to_contigs}
# Those outputs are used only by fetch_reads.
python -m falcon_unzip.mains.fetch_reads --p-ctg={input.p_ctg} --fofn={input.fofn} --ctg-list={output.ctg_list_file}
touch {output.job_done}
"""
# TODO: Proper scattering (currently in fetch_reads.py)

# This will run in 3-unzip/0-phasing/(ctg_id)/
TASK_PHASING_RUN_SCRIPT = """\

threads_aln=24

# Alignment.
aln_dir=alignment
ctg_aln_out=${{aln_dir}}/'{params.ctg_id}_sorted.bam'
mkdir -p ${{aln_dir}}
time minimap2 -a -x map-pb -t ${{threads_aln}} {input.ref_fasta} {input.read_fasta} > tmp_aln.sam

# Remove the secondary (0x100) and supplementary (0x800) alignments.
time samtools view -bS -F 0x900 tmp_aln.sam > tmp_aln.bam

time samtools sort tmp_aln.bam -o ${{ctg_aln_out}}
time samtools index ${{ctg_aln_out}}
rm tmp_aln.bam
rm tmp_aln.sam

bam_fn=${{ctg_aln_out}}
fasta_fn={input.ref_fasta}

# MAKE_HET_CALL
vmap_fn='het_call/variant_map'
vpos_fn='het_call/variant_pos'
q_id_map_fn='het_call/q_id_map.msgpack'
mkdir -p het_call
python -m falcon_unzip.mains.phasing_make_het_call --bam ${{bam_fn}} --fasta ${{fasta_fn}} --ctg-id {params.ctg_id} --vmap=${{vmap_fn}} --vpos=${{vpos_fn}} --q-id-map=${{q_id_map_fn}}

# GENERATE ASSOCIATION TABLE
atable_fn='g_atable/atable'
mkdir -p g_atable
python -m falcon_unzip.mains.phasing_generate_association_table --ctg-id {params.ctg_id} --vmap=${{vmap_fn}} --atable=${{atable_fn}}

# GET PHASED BLOCKS
phased_variant_fn='get_phased_blocks/phased_variants'
mkdir -p get_phased_blocks
python -m falcon_unzip.mains.phasing_get_phased_blocks --vmap=${{vmap_fn}} --atable=${{atable_fn}} --p-variant=${{phased_variant_fn}}

# GET PHASED READS
phased_reads_fn='get_phased_reads/phased_reads'
mkdir -p get_phased_reads
python -m falcon_unzip.mains.phasing_get_phased_reads --ctg-id={params.ctg_id} --vmap=${{vmap_fn}} --p-variant=${{phased_variant_fn}} --q-id-map=${{q_id_map_fn}} --phased-reads=${{phased_reads_fn}}

# PHASING READMAP
# TODO: read-map-dir/* as inputs
python -m falcon_unzip.mains.phasing_readmap --the-ctg-id={params.ctg_id} --rawread-ids-fn={input.rawread_ids} --pread-ids-fn={input.pread_ids} --pread-to-contigs={input.pread_to_contigs} --phased-reads=${{phased_reads_fn}} >| {output.rid_to_phase_out}.tmp
mv {output.rid_to_phase_out}.tmp {output.rid_to_phase_out}.true

mkdir -p proto
preads_ovl_dir="{params.base_dir}/1-preads_ovl"
falcon_asm_dir="{params.base_dir}/2-asm-falcon"
unzip_dir="{params.base_dir}/3-unzip"

python -m falcon_unzip.proto.extract_phased_preads \
    --ctg-id {params.ctg_id} \
    --preads ${{preads_ovl_dir}}/db2falcon/preads4falcon.fasta \
    --rid-phase-map {output.rid_to_phase_out}.true \
    --out proto/preads.fasta

ln -sf {input.ref_fasta} proto/ref.fa
time minimap2 -a -x map-pb -t ${{threads_aln}} proto/ref.fa proto/preads.fasta > proto/preads.sam

python -m falcon_unzip.proto.main_augment_pb \
    --wd ./proto/ \
    --ctg-id {params.ctg_id} \
    --p-ctg ${{falcon_asm_dir}}/p_ctg.fa \
    --p-ctg-tiling-path ${{falcon_asm_dir}}/p_ctg_tiling_path \
    --a-ctg ${{falcon_asm_dir}}/a_ctg.fa \
    --a-ctg-tiling-path ${{falcon_asm_dir}}/a_ctg_tiling_path \
    --p-variant-fn get_phased_blocks/phased_variants \
    --preads-sam proto/preads.sam \
    --extracted-ctg-fasta ${{unzip_dir}}/reads/{params.ctg_id}/ref.fa \
    --rawread-bam ${{ctg_aln_out}} \
    --rid-phase-map {output.rid_to_phase_out}.true \
    --out-updated-rid-phase_map {output.rid_to_phase_out}.tmp
mv {output.rid_to_phase_out}.tmp {output.rid_to_phase_out}
"""

TASK_PHASING_SPLIT_SCRIPT = """\
python -m falcon_unzip.mains.phasing_split --base-dir={params.topdir} --ctg-list-fn={input.ctg_list} --rawread-ids-fn={input.rawread_ids} --pread-ids-fn={input.pread_ids} --pread-to-contigs-fn={input.pread_to_contigs} --split-fn={output.split} --bash-template-fn={output.bash_template}
"""

TASK_PHASING_GATHER_SCRIPT = """\
python -m falcon_unzip.mains.phasing_gather --gathered={input.gathered} --rid-to-phase-all={output.rid_to_phase_all}
"""

TASK_HASM_SCRIPT = """\
# TODO: Needs preads.db

rm -f ./ctg_paths
python -m falcon_unzip.mains.ovlp_filter_with_phase_strict --fofn {input.las_fofn} --max-diff 120 --max-cov 120 --min-cov 1 --n-core 48 --min-len 2500 --db {input.preads_db} --rid-phase-map {input.rid_to_phase_all} > preads.p_ovl
python -m falcon_unzip.mains.phased_ovlp_to_graph preads.p_ovl --min-len 2500 > fc.log

if [[ ! -e ./ctg_paths ]]; then
    exit 1
fi

rm -f preads.p_ovl

# Create haplotigs in a safe manner.

ln -sf {input.preads4falcon} .

rm -f {output.p_ctg}

# Given sg_edges_list, utg_data, ctg_paths, preads4falcon.fasta,
# write p_ctg.fa and a_ctg_all.fa,
# plus p_ctg_tiling_path, a_ctg_tiling_path:
time python -m falcon_kit.mains.graph_to_contig

if [[ ! -e {output.p_ctg} ]]; then
    exit 1
fi
"""

TASK_GRAPH_TO_H_TIGS_SPLIT_SCRIPT = """\
asm_dir=$(dirname {input.falcon_asm_done})
hasm_dir=$(dirname {input.p_ctg})

# TODO: Should we look at ../reads/ctg_list ?

python -m falcon_unzip.mains.graphs_to_h_tigs_2 split --gathered-rid-to-phase={input.gathered_rid_to_phase} --base-dir={params.topdir} --fc-asm-path ${{asm_dir}} --fc-hasm-path ${{hasm_dir}} --rid-phase-map {input.rid_to_phase_all} --fasta {input.preads4falcon} --split-fn={output.split} --bash-template-fn={output.bash_template}

# The bash-template is just a dummy, for now.
"""

TASK_GRAPH_TO_H_TIGS_SCRIPT = """\
asm_dir=$(dirname {input.falcon_asm_done})
hasm_dir=$(dirname {input.p_ctg})

python -m falcon_unzip.mains.graphs_to_h_tigs_2 --gathered-rid-to-phase={input.gathered_rid_to_phase} --base-dir={params.topdir} --fc-asm-path ${{asm_dir}} --fc-hasm-path ${{hasm_dir}} --ctg-id all --rid-phase-map {input.rid_to_phase_all} --fasta {input.preads4falcon}

# more script -- a little bit hacky here, we should improve

#WD=$PWD
# for f in `cat ../reads/ctg_list `; do mkdir -p $WD/$f; cd $WD/$f; python -m falcon_unzip.mains.dedup_h_tigs $f; done

for f in `cat ../reads/ctg_list `
do
    mkdir -p ./$f;
    if [ -s ./$f/h_ctg.$f.fa ]
    then
        grep ">" ./$f/h_ctg.$f.fa | sed "s/^>//" >| ./$f/h_ctg_ids.$f
    else
        rm -rf ./$f/h_ctg_ids.$f
        touch ./$f/h_ctg_ids.$f
    fi
done

touch {output.htigs_done}
"""

TASK_HASM_COLLECT_SCRIPT = """\
## prepare for quviering the haplotig
## (assume we are in 3-unzip/somewhere/)

# TODO: Stop using job_done.
python -m falcon_unzip.mains.graphs_to_h_tigs_2 combine --results-fn={input.results} --done-fn={output.job_done}

find ./0-phasing -name "phased_reads" | sort | xargs cat >| all_phased_reads
#find ./2-htigs -name "h_ctg_ids.*" | sort | xargs cat >| all_h_ctg_ids
#find ./2-htigs -name "p_ctg_edges.*" | sort | xargs cat >| all_p_ctg_edges
#find ./2-htigs -name "h_ctg_edges.*" | sort | xargs cat >| all_h_ctg_edges
#find ./2-htigs -name "p_ctg.*.fa" | sort | xargs cat >| all_p_ctg.fa
#find ./2-htigs -name "h_ctg.*.fa" | sort | xargs cat >| all_h_ctg.fa

if [[ ! -s all_p_ctg.fa ]]; then
    echo "Empty all_p_ctg.fa -- No point in continuing!"
    exit 1
fi

# # Generate a GFA for only primary contigs and haplotigs.
# time python -m falcon_unzip.mains.unzip_gen_gfa_v1 --unzip-root . --p-ctg-fasta ./all_p_ctg.fa --h-ctg-fasta ./all_h_ctg.fa --preads-fasta {input.preads4falcon} >| ./asm.gfa

# # Generate a GFA of all assembly graph edges. This GFA can contain
# # edges and nodes which are not part of primary contigs and haplotigs
# time python -m falcon_unzip.mains.unzip_gen_gfa_v1 --unzip-root . --p-ctg-fasta ./all_p_ctg.fa --h-ctg-fasta ./all_h_ctg.fa --preads-fasta {input.preads4falcon} --add-string-graph >| ./sg.gfa
"""

TASK_TRACK_READS_H_SCRIPT = """\
python -m falcon_unzip.mains.get_read_hctg_map --base-dir={params.topdir} --output=./read_to_contig_map

fc_rr_hctg_track.py --stream  --db={input.r_db} --las-fofn={input.r_las_fofn} --phased-reads={input.all_phased_reads} --rawread-ids={input.rawread_ids} --read-to-contig-map=./read_to_contig_map --partials-fn=./partials.json
# That writes a msgpack partial for each raw_reads.*.las file.
# Default n_core is n_cpu/2. TODO: Configure --n-core
# (We also use a proc for LA4Falcon, so this is half.)
# It is actually limited by number of files, but in theory we could use whole machine if we have enough blocks.

abs_rawread_to_contigs=$(readlink -f {output.rawread_to_contigs}) #TODO: No readlink
cwd=$(pwd)
cd {params.topdir}
fc_rr_hctg_track2.exe --read-to-contig-map=${{cwd}}/read_to_contig_map --output=${{abs_rawread_to_contigs}} --partials-fn=${{cwd}}/partials.json
cd -

# Clean up.
python -m falcon_unzip.mains.remove_all --fofn=./partials.json

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

VC_IGNORE_ERROR='{params.vc_ignore_error}'
echo "VC_IGNORE_ERROR='$VC_IGNORE_ERROR'"
USE_BLASR='{params.use_blasr}'
echo "USE_BLASR='$USE_BLASR'"

nproc={params.pypeflow_nproc}
samtools faidx {input.ref_fasta}
if [[ $USE_BLASR != 1 ]]; then
  pbmm2 align --sort \
             {input.ref_fasta} {input.read_bam} aln-{params.ctg_id}.bam
  pbindex aln-{params.ctg_id}.bam
else
  pbalign --tmpDir=$(pwd)/tmp --nproc=$nproc --minAccuracy=0.75 --minLength=50\
          --minAnchorSize=12 --maxDivergence=30 --concordant --algorithm=blasr\
          --algorithmOptions=--useQuality --maxHits=1 --hitPolicy=random --seed=1\
            {input.read_bam} {input.ref_fasta} aln-{params.ctg_id}.bam
fi

#python -c 'import ConsensusCore2 as cc2; print cc2' # So quiver likely works.

set +e
variantCaller --algorithm=arrow -x 5 -X 120 -q 20 -j $nproc -r {input.ref_fasta} aln-{params.ctg_id}.bam\
            -o {output.cns_fasta} -o {output.cns_fastq} --minConfidence 0 -o {output.cns_vcf}
rc=$?
if [[ $rc != 0 ]]; then
    if [[ $VC_IGNORE_ERROR != 1 ]]; then
        echo ERROR variantCaller failed. Maybe no reads for this block?
        exit 1
    else
        echo WARNING variantCaller failed. Maybe no reads for this block.
        # We expect variantCaller to write files even on error, so we do not need to "touch" them.
    fi
fi
set -e

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
python -m falcon_unzip.mains.quiver_split --unzip-config-fn={input.unzip_config_fn} --p-ctg-fasta-fn={input.p_ctg_fa} --h-ctg-fasta-fn={input.h_ctg_fa} --ctg2bamfn-fn={input.ctg2bamfn} --split-fn={output.split} --bash-template-fn={output.bash_template}
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


def create_tasks_read_to_contig_map(wf, rule_writer, falcon_asm_done, raw_reads_db, preads_db, rawread_ids_fn, pread_ids_fn, read_to_contig_map_file, parameters):

    wf.addTask(gen_task(
        script=pype_tasks.TASK_DUMP_RAWREAD_IDS_SCRIPT,
        inputs={'rawread_db': raw_reads_db,
                'falcon_asm_done': falcon_asm_done,
        },
        outputs={'rawread_id_file': rawread_ids_fn,
        },
        parameters=parameters,
        rule_writer=rule_writer,
        dist=Dist(local=True), # TODO: Is this ok to run locally?
    ))

    wf.addTask(gen_task(
        script=pype_tasks.TASK_DUMP_PREAD_IDS_SCRIPT,
        inputs={'pread_db': preads_db,
                'falcon_asm_done': falcon_asm_done,
        },
        outputs={'pread_id_file': pread_ids_fn,
        },
        parameters=parameters,
        rule_writer=rule_writer,
        dist=Dist(local=True), # TODO: Is this ok to run locally?
    ))

    # From Falcon. (We will produce phased versions later.)
    sg_edges_list = '2-asm-falcon/sg_edges_list'
    utg_data = '2-asm-falcon/utg_data'
    ctg_paths = '2-asm-falcon/ctg_paths'

    inputs = {'rawread_id_file': rawread_ids_fn,
              'pread_id_file': pread_ids_fn,
              'sg_edges_list': sg_edges_list,
              'utg_data': utg_data,
              'ctg_paths': ctg_paths}
    wf.addTask(gen_task(
        script=pype_tasks.TASK_GENERATE_READ_TO_CTG_MAP_SCRIPT,
        inputs=inputs,
        outputs={'read_to_contig_map': read_to_contig_map_file},
        parameters=parameters,
        rule_writer=rule_writer,
        dist=Dist(local=True), # TODO: Is this ok to run locally?
    ))


def run_workflow(wf, config, unzip_config_fn, rule_writer):
    default_njobs = int(config['job.defaults']['njobs'])
    wf.max_jobs = default_njobs

    falcon_asm_done_fn = './2-asm-falcon/falcon_asm_done'
    p_ctg_fn = './2-asm-falcon/p_ctg.fa'
    raw_reads_db_fn = './0-rawreads/build/raw_reads.db'
    preads_db_fn = './1-preads_ovl/build/preads.db'
    r_las_fofn_fn = './0-rawreads/las-merge-combine/las_fofn.json'
    p_las_fofn_fn = './1-preads_ovl/las-merge-combine/las_fofn.json'

    read_to_contig_map_fn = '3-unzip/reads/get_read_ctg_map/read_to_contig_map'
    rawread_ids_fn = '3-unzip/reads/dump_rawread_ids/rawread_ids'
    pread_ids_fn = '3-unzip/reads/dump_pread_ids/pread_ids'
    # This has lots of inputs from falcon stages 0, 1, and 2.
    create_tasks_read_to_contig_map(wf, rule_writer, falcon_asm_done_fn, raw_reads_db_fn, preads_db_fn, rawread_ids_fn, pread_ids_fn, read_to_contig_map_fn, {})

    ctg_list_fn = './3-unzip/reads/ctg_list'
    rawread_to_contigs_fn = './3-unzip/reads/rawread_to_contigs'
    pread_to_contigs_fn = './3-unzip/reads/pread_to_contigs'

    LOG.info('config=\n {}'.format(config))
    Unzip_config = config['Unzip']
    fasta_fofn_fn = Unzip_config.get('input_fofn') #, './input.fofn') # from user config, usually

    wf.addTask(gen_task(
            script=TASK_TRACK_READS_SCRIPT,
            inputs={
                'fofn': fasta_fofn_fn,
                'read_to_contig_map': read_to_contig_map_fn,
                'p_ctg': p_ctg_fn,
                'r_db': raw_reads_db_fn,
                'p_db': preads_db_fn,
                'r_las_fofn': r_las_fofn_fn,
                'p_las_fofn': p_las_fofn_fn,
            },
            outputs={
                'job_done': './3-unzip/reads/track_reads_done',
                'ctg_list_file': ctg_list_fn,
                'rawread_to_contigs': rawread_to_contigs_fn,
                'pread_to_contigs': pread_to_contigs_fn,
            },
            parameters={},
            rule_writer=rule_writer,
            dist=Dist(
                #NPROC=4, # no longer hard-coded
                job_dict=config['job.step.unzip.track_reads'],
                use_tmpdir=False,
            ),
    ))

    phasing_all_units_fn = './3-unzip/0-phasing/phasing-split/all-units-of-work.json'
    phasing_run_bash_template_fn ='./3-unzip/0-phasing/phasing-split/bash-template.sh'

    wf.addTask(gen_task(
        script=TASK_PHASING_SPLIT_SCRIPT,
        inputs=dict(
            ctg_list=ctg_list_fn,
            rawread_ids=rawread_ids_fn,
            pread_ids=pread_ids_fn,
            pread_to_contigs=pread_to_contigs_fn,
        ),
        outputs=dict(
            split=phasing_all_units_fn,
            bash_template=phasing_run_bash_template_fn,
        ),
        parameters={},
        rule_writer=rule_writer,
        dist=Dist(local=True),
    ))

    gathered_rid_to_phase_fn = './3-unzip/0-phasing/gathered-rid-to-phase/gathered.json'

    gen_parallel_tasks(
        wf, rule_writer,
        phasing_all_units_fn, gathered_rid_to_phase_fn,
        run_dict=dict(
            bash_template_fn=phasing_run_bash_template_fn,
            script=TASK_PHASING_RUN_SCRIPT,
            inputs={
                'units_of_work': './3-unzip/0-phasing/phasing-chunks/{ctg_id}/some-units-of-work.json',
            },
            outputs={
                'results': './3-unzip/0-phasing/{ctg_id}/phasing-result-list.json',
            },
            parameters={},
        ),
        dist=Dist(NPROC=24, job_dict=config['job.step.unzip.blasr_aln']),
    )

    concatenated_rid_to_phase_fn = './3-unzip/1-hasm/concatenated-rid-to-phase/rid_to_phase.all'

    wf.addTask(gen_task(
        script=TASK_PHASING_GATHER_SCRIPT,
        inputs={'gathered': gathered_rid_to_phase_fn,
        },
        outputs={'rid_to_phase_all': concatenated_rid_to_phase_fn,
        },
        parameters={},
        rule_writer=rule_writer,
        dist=Dist(local=True),
    ))

    preads4falcon_fn = './1-preads_ovl/db2falcon/preads4falcon.fasta'

    hasm_p_ctg_fn = './3-unzip/1-hasm/p_ctg.fa'
    wf.addTask(gen_task(
            script=TASK_HASM_SCRIPT,
            inputs={
                'preads_db': preads_db_fn,
                'preads4falcon': preads4falcon_fn,
                'las_fofn': p_las_fofn_fn,
                'rid_to_phase_all': concatenated_rid_to_phase_fn,
            },
            outputs={
                'p_ctg': hasm_p_ctg_fn,
            },
            parameters={},
            rule_writer=rule_writer,
            dist=Dist(NPROC=1, job_dict=config['job.step.unzip.hasm']),
    ))

    #htigs_done_fn = './3-unzip/2-htigs/htigs.done'
    g2h_all_units_fn = './3-unzip/2-htigs/split/all-units-of-work.json'
    dummy_fn = './3-unzip/2-htigs/split/dummy.sh'
    wf.addTask(gen_task(
            script=TASK_GRAPH_TO_H_TIGS_SPLIT_SCRIPT,
            inputs={
                'falcon_asm_done': falcon_asm_done_fn,
                'preads4falcon': preads4falcon_fn,
                'rid_to_phase_all': concatenated_rid_to_phase_fn,
                'gathered_rid_to_phase': gathered_rid_to_phase_fn,
                'p_ctg': hasm_p_ctg_fn,
            },
            outputs={
                'split': g2h_all_units_fn,
                'bash_template': dummy_fn,
            },
            parameters={},
            rule_writer=rule_writer,
            dist=Dist(NPROC=1, job_dict=config['job.step.unzip.hasm']),
    ))
    TASK_GTOH_APPLY_UNITS_OF_WORK = """\
    python -m falcon_unzip.mains.graphs_to_h_tigs_2 apply --units-of-work-fn={input.units_of_work} --results-fn={output.results}

    #--bash-template-fn= # not needed
    """
    gathered_g2h_fn = './3-unzip/2-htigs/gathered/gathered.json'
    gen_parallel_tasks(
        wf, rule_writer,
        g2h_all_units_fn, gathered_g2h_fn,
        run_dict=dict(
            bash_template_fn=dummy_fn,
            script='DUMMY',
            inputs={
                'units_of_work': './3-unzip/2-htigs/chunks/{chunk_id}/some-units-of-work.json',
            },
            outputs={
                'results': './3-unzip/2-htigs/{chunk_id}/result-list.json',
            },
            parameters={},
        ),
        dist=Dist(
            NPROC=24,
            job_dict=config['job.step.unzip.blasr_aln'],
            use_tmpdir=False,
        ),
        run_script=TASK_GTOH_APPLY_UNITS_OF_WORK,
    )

    job_done = './3-unzip/hasm_done'
    wf.addTask(gen_task(
            script=TASK_HASM_COLLECT_SCRIPT,
            inputs={
                'preads4falcon': preads4falcon_fn,
                'results': gathered_g2h_fn,
            },
            outputs={
                'job_done': job_done,
                'all_phased_reads': './3-unzip/all_phased_reads',
                'p_ctg_fa': './3-unzip/all_p_ctg.fa',
                'h_ctg_fa': './3-unzip/all_h_ctg.fa',
            },
            parameters={},
            rule_writer=rule_writer,
            dist=Dist(
                NPROC=1,
                job_dict=config['job.step.unzip.hasm'],
                use_tmpdir=False,
            ),
    ))

    ################
    # 4-polish stage

    input_bam_fofn = os.path.relpath(config['Unzip']['input_bam_fofn']) # All input paths should be relative, for snakemake.
    rawreads_db_fn = './0-rawreads/build/raw_reads.db'
    #preads_db_fn = './1-preads_ovl/build/preads.db'
    r_las_fofn_fn = './0-rawreads/las-merge-combine/las_fofn.json'
    #p_las_fofn_fn = './1-preads_ovl/las-merge-combine/las_fofn.json'
    all_phased_reads_fn = './3-unzip/all_phased_reads'
    rawread_ids_fn = './3-unzip/reads/dump_rawread_ids/rawread_ids'

    track_reads_rr2c = './4-polish/track-reads/rawread_to_contigs'

    wf.addTask(gen_task(
        script=TASK_TRACK_READS_H_SCRIPT,
        inputs={
            #'input_bam_fofn': input_bam_fofn,
            'r_db': rawreads_db_fn,
            'r_las_fofn': r_las_fofn_fn,
            'all_phased_reads': all_phased_reads_fn,
            'rawread_ids': rawread_ids_fn,
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

    read2ctg = './4-polish/select-reads/read2ctg.msgpack'
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

    merged_fofn = './4-polish/merge-reads/merged.fofn'
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

    segr_all_units_fn ='./4-polish/segregate-split/all-units-of-work.json'
    segr_run_bash_template_fn ='./4-polish/segregate-split/bash-template.sh'
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
    gathered_fn = './4-polish/segregate-gathered/segregated-bam.json'
    gen_parallel_tasks(
        wf, rule_writer,
        segr_all_units_fn, gathered_fn,
        run_dict=dict(
            bash_template_fn=segr_run_bash_template_fn,
            script=TASK_SEGREGATE_RUN_SCRIPT,
            inputs={
                'units_of_work': './4-polish/segregate-chunks/{segr}/some-units-of-work.json',
            },
            outputs={
                'results': './4-polish/segregate-run/{segr}/segregated-bam-fn-list.json',
            },
            parameters={},
        ),
        dist=Dist(NPROC=4, MB=16, # actually NPROC=1, but our qsub jobs rarely report mem needs
            job_dict=config['job.step.unzip.track_reads'],
            use_tmpdir=False,
        ),
    )

    ctg2segregated_bamfn = './4-polish/segregated-bam/ctg2segregated_bamfn.msgpack'
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

    quiver_all_units_fn ='./4-polish/quiver-split/all-units-of-work.json'
    quiver_run_bash_template_fn ='./4-polish/quiver-split/bash-template.sh'
    wf.addTask(gen_task(
        script=TASK_QUIVER_SPLIT_SCRIPT,
        inputs={
                'p_ctg_fa': './3-unzip/all_p_ctg.fa',
                'h_ctg_fa': './3-unzip/all_h_ctg.fa',
                'ctg2bamfn': ctg2segregated_bamfn,
                'unzip_config_fn': unzip_config_fn,
        },
        outputs={
                'split': quiver_all_units_fn,
                'bash_template': quiver_run_bash_template_fn,
        },
        parameters={},
        rule_writer=rule_writer,
        #dist=Dist(local=True),
        # lots of fasta parsing, so must not run locally
        dist=Dist(NPROC=1,
            job_dict=config['job.step.unzip.quiver'],
        ),
    ))

    int_gathered_fn = '4-polish/cns-gather/intermediate/int.gathered.json'
    gen_parallel_tasks(
        wf, rule_writer,
        quiver_all_units_fn, int_gathered_fn,
        run_dict=dict(
            bash_template_fn=quiver_run_bash_template_fn,
            script='',
            inputs={
                'units_of_work': './4-polish/quiver-chunks/{ctg_id}/some-units-of-work.json',
            },
            outputs={
                'results': './4-polish/quiver-run/{ctg_id}/results.json',
            },
        ),
        dist=Dist(
            job_dict=config['job.step.unzip.quiver'],
        ),
    )
    gathered_quiver = '4-polish/cns-gather/gathered_quiver.json'
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
            'cns_p_ctg_fasta': '4-polish/cns-output/cns_p_ctg.fasta',
            'cns_p_ctg_fastq': '4-polish/cns-output/cns_p_ctg.fastq',
            'cns_h_ctg_fasta': '4-polish/cns-output/cns_h_ctg.fasta',
            'cns_h_ctg_fastq': '4-polish/cns-output/cns_h_ctg.fastq',
            'job_done': '4-polish/cns-output/job_done',
        },
        rule_writer=rule_writer,
        dist=Dist(NPROC=1),
    ))

    wf.refreshTargets()
