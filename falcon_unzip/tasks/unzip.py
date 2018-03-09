from pypeflow.simple_pwatcher_bridge import (
    makePypeLocalFile, fn,
    PypeTask,
)
from .pype import (gen_task, gen_parallel_tasks)
from .. import io
from falcon_kit import pype_tasks
import logging
import os

LOG = logging.getLogger(__name__)


TASK_TRACK_READS_SCRIPT = """\
# Also require read_to_contig_map.
python -m falcon_unzip.mains.rr_ctg_track --base-dir={params.topdir} --output=rawread_to_contigs
python -m falcon_unzip.mains.pr_ctg_track --base-dir={params.topdir} --output=pread_to_contigs
# Those outputs are used only by fetch_reads.
python -m falcon_unzip.mains.fetch_reads --base-dir={params.topdir} --fofn={input.fofn} --ctg-list={output.ctg_list_file}
touch {output.job_done}
"""
# TODO: Proper scattering (currently in fetch_reads.py)

# This will run in 3-unzip/0-phasing/(ctg_id)/
TASK_RUN_PHASING_SCRIPT = """\

# BLASR
ctg_aln_out='blasr/{params.ctg_id}_sorted.bam'
mkdir -p blasr
time blasr {input.read_fasta} {input.ref_fasta} --noSplitSubreads --clipping subread\
 --hitPolicy randombest --randomSeed 42 --bestn 1 --minPctIdentity 70.0\
 --minMatch 12  --nproc 24 --bam --out tmp_aln.bam
#samtools view -bS tmp_aln.sam | samtools sort - {params.ctg_id}_sorted
samtools sort tmp_aln.bam -o ${{ctg_aln_out}}
samtools index ${{ctg_aln_out}}
rm tmp_aln.bam

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
python -m falcon_unzip.mains.phasing_readmap --the-ctg-id={params.ctg_id} --read-map-dir=../../reads --phased-reads=${{phased_reads_fn}} >| {output.rid_to_phase_out}.tmp
mv {output.rid_to_phase_out}.tmp {output.rid_to_phase_out}
"""

TASK_PHASING_SCATTER_SCRIPT = """\
python -m falcon_unzip.mains.phasing_scatter --ctg-list-fn={input.ctg_list} --scattered-fn={output.scattered}
"""

TASK_GET_RID_TO_PHASE_ALL_SCRIPT = """\
rm -f {output.rid_to_phase_all}
for fn in {input}; do
  cat $fn >> {output.rid_to_phase_all}
done
"""

TASK_PHASING_GATHER_SCRIPT = """\
python -m falcon_unzip.mains.phasing_gather --gathered={input.gathered} --rid-to-phase-all={output.rid_to_phase_all}
"""


TASK_HASM_SCRIPT = """\
python -m falcon_unzip.mains.ovlp_filter_with_phase --fofn {input.las_fofn} --max_diff 120 --max_cov 120 --min_cov 1 --n_core 48 --min_len 2500 --db ../../1-preads_ovl/preads.db --rid_phase_map {input.rid_to_phase_all} > preads.p_ovl
python -m falcon_unzip.mains.phased_ovlp_to_graph preads.p_ovl --min_len 2500 > fc.log
if [ -e ../../1-preads_ovl/preads4falcon.fasta ];
then
  ln -sf ../../1-preads_ovl/preads4falcon.fasta .
else
  ln -sf ../../1-preads_ovl/db2falcon/preads4falcon.fasta .
fi
python -m falcon_unzip.mains.graphs_to_h_tigs --fc_asm_path ../../2-asm-falcon/ --fc_hasm_path ./ --ctg_id all --rid_phase_map {input.rid_to_phase_all} --fasta preads4falcon.fasta

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

cd $WD
touch {output.job_done}
"""


def create_tasks_read_to_contig_map(wf, rule_writer, read_to_contig_map_file):
    falcon_asm_done = './2-asm-falcon/falcon_asm_done'

    rawread_db = '0-rawreads/raw_reads.db'
    rawread_ids = '3-unzip/reads/dump_rawread_ids/rawread_ids'

    wf.addTask(gen_task(
        script=pype_tasks.TASK_DUMP_RAWREAD_IDS_SCRIPT,
        inputs={'rawread_db': rawread_db,
                'falcon_asm_done': falcon_asm_done,
        },
        outputs={'rawread_id_file': rawread_ids,
        },
        parameters={},
        rule_writer=rule_writer,
    ))

    pread_db = '1-preads_ovl/preads.db'
    pread_ids = '3-unzip/reads/dump_pread_ids/pread_ids'

    wf.addTask(gen_task(
        script=pype_tasks.TASK_DUMP_PREAD_IDS_SCRIPT,
        inputs={'pread_db': pread_db,
                'falcon_asm_done': falcon_asm_done,
        },
        outputs={'pread_id_file': pread_ids,
        },
        parameters={},
        rule_writer=rule_writer,
    ))

    sg_edges_list = '2-asm-falcon/sg_edges_list'
    utg_data = '2-asm-falcon/utg_data'
    ctg_paths = '2-asm-falcon/ctg_paths'

    inputs = {'rawread_id_file': rawread_ids,
              'pread_id_file': pread_ids,
              'sg_edges_list': sg_edges_list,
              'utg_data': utg_data,
              'ctg_paths': ctg_paths}
    wf.addTask(gen_task(
        script=pype_tasks.TASK_GENERATE_READ_TO_CTG_MAP_SCRIPT,
        inputs=inputs,
        outputs={'read_to_contig_map': read_to_contig_map_file},
        parameters={},
        rule_writer=rule_writer,
    ))


def run_workflow(wf, config, rule_writer):
    parameters = {
        'sge_option': config['sge_track_reads'],
        #'topdir': os.getcwd(),
    }
    read_to_contig_map_file = '3-unzip/reads/get_read_ctg_map/read_to_contig_map'
    # This has lots of inputs from falcon stages 0, 1, and 2.
    create_tasks_read_to_contig_map(wf, rule_writer, read_to_contig_map_file)

    ctg_list_file = './3-unzip/reads/ctg_list'
    fofn_file = config.get('input_fofn', './input.fofn') # from user config, usually

    wf.addTask(gen_task(
            script=TASK_TRACK_READS_SCRIPT,
            inputs={
                'fofn': fofn_file,
                'read_to_contig_map': read_to_contig_map_file,
            },
            outputs={
                'job_done': './3-unzip/reads/track_reads_done',
                'ctg_list_file': ctg_list_file,
            },
            parameters=parameters,
            rule_writer=rule_writer,
    ))

    scattered = './3-unzip/1-hasm/scattered/scattered.json'

    wf.addTask(gen_task(
        script=TASK_PHASING_SCATTER_SCRIPT,
        inputs=dict(
            ctg_list=ctg_list_file,
        ),
        outputs=dict(
            scattered=scattered,
        ),
        parameters={},
        rule_writer=rule_writer,
    ))

    parameters = {#'job_uid': 'aln-' + ctg_id, 'ctg_id': ctg_id,
                    'sge_option': config['sge_blasr_aln'],
    }
    gathered = './3-unzip/1-hasm/gathered-rid-to-phase/gathered.json'

    gen_parallel_tasks(
        wf, rule_writer,
        scattered, gathered,
        run_dict=dict(
            script=TASK_RUN_PHASING_SCRIPT,
            inputs={
                'ref_fasta': './3-unzip/reads/{ctg_id}/ref.fa',
                'read_fasta': './3-unzip/reads/{ctg_id}/reads.fa',
            },
            outputs={
                'rid_to_phase_out': './3-unzip/0-phasing/{ctg_id}/rid_to_phase',
            },
            parameters=parameters,
        ),
    )

    concatenated_rid_to_phase_file = './3-unzip/1-hasm/concatenated-rid-to-phase/rid_to_phase.all'

    wf.addTask(gen_task(
        script=TASK_PHASING_GATHER_SCRIPT,
        inputs={'gathered': gathered,
        },
        outputs={'rid_to_phase_all': concatenated_rid_to_phase_file,
        },
        parameters={},
        rule_writer=rule_writer,
    ))

    parameters = {
            'sge_option': config['sge_hasm'],
    }
    las_fofn_file = './1-preads_ovl/las-gather/las_fofn.json'
    job_done = './3-unzip/1-hasm/hasm_done'

    wf.addTask(gen_task(
            script=TASK_HASM_SCRIPT,
            inputs={
                'rid_to_phase_all': concatenated_rid_to_phase_file,
                'las_fofn': las_fofn_file,
            },
            outputs={
                'job_done': job_done,
            },
            parameters=parameters,
            rule_writer=rule_writer,
    ))
    unzip_phasing_concurrent_jobs = config['unzip_phasing_concurrent_jobs']
    wf.max_jobs = unzip_phasing_concurrent_jobs

    wf.refreshTargets()
