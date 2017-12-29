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
python -m falcon_unzip.mains.fetch_reads --base-dir={params.topdir} --fofn={output.fofn}
touch {output.job_done}
# Also produce ctg_list_file.
"""


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


def task_hasm(self):
    rid_to_phase_all = fn(self.rid_to_phase_all)
    las_fofn = fn(self.las_fofn)
    job_done = fn(self.job_done)
    #config = self.parameters['config']

    script_fn = 'hasm.sh'
    script = """\
set -vex
trap 'touch {job_done}.exit' EXIT
hostname
date

python -m falcon_unzip.mains.ovlp_filter_with_phase --fofn {las_fofn} --max_diff 120 --max_cov 120 --min_cov 1 --n_core 48 --min_len 2500 --db ../../1-preads_ovl/preads.db --rid_phase_map {rid_to_phase_all} > preads.p_ovl
python -m falcon_unzip.mains.phased_ovlp_to_graph preads.p_ovl --min_len 2500 > fc.log
if [ -e ../../1-preads_ovl/preads4falcon.fasta ];
then
  ln -sf ../../1-preads_ovl/preads4falcon.fasta .
else
  ln -sf ../../1-preads_ovl/db2falcon/preads4falcon.fasta .
fi
python -m falcon_unzip.mains.graphs_to_h_tigs --fc_asm_path ../../2-asm-falcon/ --fc_hasm_path ./ --ctg_id all --rid_phase_map {rid_to_phase_all} --fasta preads4falcon.fasta

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

cd ../
date
touch {job_done}
""".format(**locals())

    with open(script_fn, 'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn


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


def get_hasm_task(config, gathered_rid_to_phase_file, las_fofn_file, job_done):
    parameters = {
            'sge_option': config['sge_hasm'],
    }
    make_hasm_task = PypeTask(
            inputs={
                'rid_to_phase_all': gathered_rid_to_phase_file,
                'las_fofn': las_fofn_file,
            },
            outputs={
                'job_done': job_done,
            },
            parameters=parameters,
            )
    return make_hasm_task(task_hasm)


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
                'job_done': './3-unzip/reads/track_reads_done', # TODO: Depend directly on ctg_list_file instead.
                'ctg_list_file': ctg_list_file,
            },
            parameters=parameters,
            rule_writer=rule_writer,
    ))

    # Refresh so that ctg_list_file is available. TODO: Proper scattering.
    wf.refreshTargets()

    gathered_rid_to_phase_file = makePypeLocalFile('./3-unzip/1-hasm/gathered-rid-to-phase/rid_to_phase.all')
    phasing_tasks = list(create_phasing_tasks(config, ctg_list_file, gathered_rid_to_phase_file))
    wf.addTasks(phasing_tasks)

    las_fofn_file = makePypeLocalFile('./1-preads_ovl/merge-gather/las.fofn') #'2-asm-falcon/las.fofn'
    job_done = makePypeLocalFile('./3-unzip/1-hasm/hasm_done')
    hasm_task = get_hasm_task(config, gathered_rid_to_phase_file, las_fofn_file, job_done)
    wf.addTask(hasm_task)

    unzip_phasing_concurrent_jobs = config['unzip_phasing_concurrent_jobs']
    wf.max_jobs = unzip_phasing_concurrent_jobs
    wf.refreshTargets()
