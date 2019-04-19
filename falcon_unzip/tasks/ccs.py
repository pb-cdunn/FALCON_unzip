from __future__ import absolute_import
from falcon_kit.pype import Dist
#from falcon_kit import pype_tasks
from falcon_kit.pype import (wrap_gen_task as gen_task, gen_parallel_tasks, Dist)
from .. import io
import glob
import logging
import os

LOG = logging.getLogger(__name__)

TASK_POLISH="""
## gathering phased and unphased reads
samtools faidx {input.FA} {params.ctg} > ref.fasta
perl -lane 'print $F[0] if $F[1] eq "{params.ctg}" && $F[2] != -1' {input.READTOCTG} | sort | uniq > readnames.txt
samtools view -F 1796 {input.UBAM} {params.ctg} | cut -f 1 | sort | uniq >> readnames.txt
samtools fqidx -r readnames.txt {input.FQ} > reads.fastq
minimap2 -a -x asm5 -t 1 ref.fasta reads.fastq | samtools view -F 1796 > aln.sam
racon reads.fastq aln.sam ref.fasta > {output.POL}
"""

TASK_GATHER_UNPHASED="""
find {input.ALL} > bam_list.txt
samtools merge -b bam_list.txt tmp.bam
samtools sort -o {output.UBAM} tmp.bam
samtools index {output.UBAM}
"""

TASK_PLACE_UNPHASED="""
python -m falcon_unzip.mains.polish_unphased_readmapping --ctg {params.ctg} --fai {input.fai} --out-read-names read_names.txt --out-ref-names ref_names.txt --read-to-ctg {input.readtoctg}
samtools fqidx -r read_names.txt {input.fq} > reads.fastq
samtools faidx -r ref_names.txt {input.fa}  > ref.fasta
minimap2 -a -x asm5 -t 1 ref.fasta reads.fastq | samtools view -bS -F 1796 > aln.bam
"""

TASK_POLISH_PREAMBLE="""
cat {input.P} {input.H} > {output.COMBINED}
samtools faidx {output.COMBINED}
cat ../../3-unzip/2-htigs/chunk_*/uow-*/*ctg_edges* > {output.EDGES}
python -m falcon_unzip.mains.polish_read_to_ctg --rid-to-phase-fn {input.RIDTOPHASE} --edges-fn {output.EDGES} --lookup-fn {input.READNAMELOOKUP} --out {output.READTOCTG}

# Assume FQ is external. (TODO: Use a program so we can be more lenient.)
ln -sf {input.FQ} {output.FQO}
samtools fqidx {output.FQO}
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

TASK_GRAPH_TO_H_TIGS_SPLIT_SCRIPT = """\
asm_dir=$(dirname {input.falcon_asm_done})
hasm_dir=$(dirname {input.p_ctg})

# TODO: Should we look at ../reads/ctg_list ?

python -m falcon_unzip.mains.graphs_to_h_tigs_2 split --gathered-rid-to-phase={input.gathered_rid_json} --base-dir={params.topdir} --fc-asm-path ${{asm_dir}} --fc-hasm-path ${{hasm\
_dir}} --rid-phase-map {input.rid_to_phase_all} --fasta {input.preads4falcon} --split-fn={output.split} --bash-template-fn={output.bash_template}

# The bash-template is just a dummy, for now.
"""


TASK_HASM_SCRIPT = """\
# TODO: Needs preads.db

rm -f ./ctg_paths
python -m falcon_unzip.mains.ovlp_filter_with_phase_strict --fofn {input.las_fofn} --max-diff 120 --max-cov 120 --min-cov 1 --n-core 48 --min-len 2500 --db {input.preads_db} --rid-phas\
e-map {input.rid_to_phase_all} > preads.p_ovl
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

TASK_PHASING_GATHER_SCRIPT = """\
# creates a master table of rid to phase
cat {input.ALL} > {output.rid_to_phase_all}

# creates the needed gathering JSON
find {input.ALL} | xargs -I [] readlink -f [] | python -m falcon_unzip.mains.gen_rid_gathered_json > {output.gathered_rid_json}
"""

TASK_READ_PHASING = """

        #TODO: break up command, and maybe remove some deps.

        python -m falcon_unzip.mains.phasing_make_het_call --ctg-id {params.ctg} --bam-fn {input.BAM} --fasta-fn {input.T} --vmap-fn het_calls.ctg.vmap --vpos-fn het_calls.ctg.vpos --q-id-map-fn het_calls.ctg.msgpack
        python -m falcon_unzip.mains.phasing_generate_association_table --ctg-id {params.ctg} --vmap=het_calls.ctg.vmap --atable=association_table.ctg.astab
        python -m falcon_unzip.mains.phasing_get_phased_blocks --vmap=het_calls.ctg.vmap --atable=association_table.ctg.astab --p-variant=phased_vars.ctg.phased.txt
        python -m falcon_unzip.mains.phasing_get_phased_reads --ctg-id={params.ctg} --vmap=het_calls.ctg.vmap --p-variant=phased_vars.ctg.phased.txt --q-id-map=het_calls.ctg.msgpack --phased-reads=phased_reads.ctg.phased.txt

        #reformats the data keeping the last, second..forth columns
        cat phased_reads.ctg.phased.txt | perl -lane 'print "$F[-1] $F[1] $F[2] $F[3]"' >|   phased_reads.ctg.phased.reads.reformat.txt

        mkdir -p proto

        #pulls the reference into the proto dir
        samtools faidx {input.T} {params.ctg} > proto/ref.fa

        python -m falcon_unzip.proto.main_augment_pb --wd ./proto/ --ctg-id {params.ctg}     --p-ctg {input.PCTG} --p-ctg-tiling-path {input.PTILE} --a-ctg {input.ACTG} --a-ctg-tiling-path {input.ATILE}  --p-variant-fn phased_vars.ctg.phased.txt --preads-sam {input.BAM}  --extracted-ctg-fasta {input.T} --rawread-bam {input.BAM}  --rid-phase-map phased_reads.ctg.phased.reads.reformat.txt  --out-updated-rid-phase_map rid_to_phase.tmp

        #grabs the names of all the CCS reads that associate with this primary contig.
        grep -w {params.ctg} {input.RID_TO_CTG} > rid_to_ctg.txt

        #converts the CCS read names into DAZDB read ids.
        python -m falcon_unzip.mains.db_to_ccs_id --lookup {input.readname_lookup} --rid-to-phase rid_to_phase.tmp --rid-to-ctg rid_to_ctg.txt --output {output.M} --ctg {params.ctg}

        python -m falcon_unzip.proto.extract_phased_preads --ctg-id {params.ctg} --preads ../../../../1-preads_ovl/db2falcon/preads4falcon.fasta --rid-phase-map {output.M} --out proto/preads.fasta

        time minimap2 -a -x map-pb -t 1 proto/ref.fa proto/preads.fasta > proto/preads.sam
"""


TASK_READNAME_LOOKUP = """
       #TODO: move into a script rather than a hard to grok command line.
       #This command makes a two column file where one column is the CCS read name and the other is DAZDB ID.


       paste <(DBdump -hr ../../1-preads_ovl/build/preads.db | grep "^L" ) <(DBdump -hr ../../1-preads_ovl/build/preads.db | grep "^H" ) | perl -lane '$ln = sprintf("%09d", $. -1); @Z = split /\s+/, $_; print "$ln\t$Z[-1]/$Z[1]/ccs"' > {output.readname_lookup}
"""

TASK_SORT = """
        samtools sort -@ {params.pypeflow_nproc} -o {output.OBAMS} {input.IBAMA}
        samtools index {output.OBAMS}
        samtools view -F 3844 {output.OBAMS}  | perl -lane '$F[2] =~ s/\-.*//; print "$F[0] $F[2]"' > {output.RID_TO_CTG}
"""
TASK_MAP = """
        minimap2 -a -x asm5 -t {params.pypeflow_nproc} {input.T} {input.R} | samtools view -F 3840 -bS > {output.OBAMA}
"""
TASK_PREAMBLE = """
        cat {input.P} {input.A} > {output.FA}
        samtools faidx {output.FA}
        samtools faidx {input.P}
"""



def run_workflow(wf, config, unzip_config_fn, rule_writer):
    default_njobs = int(config['job.defaults']['njobs'])
    wf.max_jobs = default_njobs

    #LOG.info('config=\n {}'.format(config))
    Unzip_config = config['Unzip']

    falcon_asm_done_fn = './2-asm-falcon/falcon_asm_done'

    #ifastq_fn = '/home/zkronenberg/dump/hg002_chr6_dataset/hg002_chr6_size_filt.fq'
    ifastq_fn = Unzip_config['fastq']
    LOG.info('Input fastq="{}"'.format(ifastq_fn))

    asm_dir = './2-asm-falcon'
    p_ctg_fn = os.path.join(asm_dir, 'p_ctg.fa')
    a_ctg_fn = os.path.join(asm_dir, 'a_ctg.fa')
    p_tile_fn = os.path.join(asm_dir, 'p_ctg_tiling_path')
    a_tile_fn = os.path.join(asm_dir, 'a_ctg_tiling_path')
    aln_cpu = '16'
    sort_cpu = '3'

    # For now, use the same job-distribution parameters everywhere.
    dist = Dist(
        job_dict=config['job.defaults'],
        #NPROC=4,
        #use_tmpdir=False,
    )

    # For strictly local jobs, use this.
    distl = Dist(
        #job_dict=config['job.defaults'],
        NPROC=1,
        local=True,
        use_tmpdir=False,
    )

    p_ctg_fai_fn = "./2-asm-falcon/p_ctg.fa.fai"

    wf.addTask(gen_task(
            script=TASK_PREAMBLE,
            inputs={
                "P": p_ctg_fn,
                "A": a_ctg_fn,
            },
            outputs={
                "FA": "3-unzip/ctgs/concat.fa",
                "FAI": "3-unzip/ctgs/concat.fa.fai",
            },
            parameters={},
            dist = Dist(
                job_dict=config['job.defaults'],
                NPROC=1,
            )
    ))

    wf.refreshTargets()

    CTGS = []

    with open(p_ctg_fai_fn) as f:
          for line in f:
            lc = line.strip().split("\t")
            if(lc[0][0] == "0"):
                CTGS.append(lc[0])


    rid_to_ctg = "./3-unzip/sorting/rid_to_cgt.txt"

    wf.addTask(gen_task(
            script=TASK_MAP,
            inputs={
                "T": "3-unzip/ctgs/concat.fa",
                "R": ifastq_fn,
            },
            outputs={
                "OBAMA": "3-unzip/mapping/reads_mapped.bam",
            },
            parameters={},
            dist = Dist(
                job_dict=config['job.defaults'],
                NPROC=aln_cpu, #TODO
            )
    ))


    wf.addTask(gen_task(
            script=TASK_SORT,
            inputs={
                "IBAMA"      : "3-unzip/mapping/reads_mapped.bam",
            },
            outputs={
                "OBAMS"  : "3-unzip/sorting/reads_mapped.sorted.bam",
                "OBAMSAI": "3-unzip/sorting/reads_mapped.sorted.bam.bai",
                "RID_TO_CTG" : rid_to_ctg,
            },
            parameters={},
            dist = Dist(
                job_dict=config['job.defaults'],
                NPROC=sort_cpu, #TODO
            )
    ))

    readname_lookup = "3-unzip/readnames/readname_lookup.txt"

    wf.addTask(gen_task(
            script=TASK_READNAME_LOOKUP,
            inputs={
                "OBAMSAI": "3-unzip/sorting/reads_mapped.sorted.bam.bai",
            },
            outputs={
                'readname_lookup'  : readname_lookup,
            },
            parameters={},
            dist = Dist(
                job_dict=config['job.defaults'],
                NPROC=sort_cpu, #TODO
            )
    ))

    collected = dict()

    for ctg in CTGS:
        rid_to_phase_fn = "3-unzip/0-phasing/{}/uow-fake/rid_to_phase".format(ctg)
        collected['ctg' + ctg] = rid_to_phase_fn
        wf.addTask(gen_task(
            script=TASK_READ_PHASING,
            inputs={
                "PCTG": p_ctg_fn,
                "ACTG": a_ctg_fn,
                "PTILE": p_tile_fn,
                "ATILE": a_tile_fn,
                "BAM": "3-unzip/sorting/reads_mapped.sorted.bam",
                "T": "3-unzip/ctgs/concat.fa",
                "readname_lookup" : readname_lookup,
                "RID_TO_CTG"      : rid_to_ctg,
            },
            outputs={
                "M": rid_to_phase_fn,
            },
            parameters={
                'ctg': ctg,
            },
            dist = Dist(
                job_dict=config['job.defaults'],
                #NPROC=???, #TODO
            )
        ))

    wf.refreshTargets()

    check(sorted(collected.itervalues()), sorted(glob.glob('3-unzip/0-phasing/*/uow-fake/rid_to_phase')))

    concatenated_rid_to_phase_fn = "3-unzip/0-phasing/gathered-rid-to-phase/rid_to_phase.all"
    gathered_rid_to_phase_json   = "3-unzip/0-phasing/gathered-rid-to-phase/gathered.json"

    wf.addTask(gen_task(
        script=TASK_PHASING_GATHER_SCRIPT,
        inputs=collected,
        outputs={'rid_to_phase_all'  : concatenated_rid_to_phase_fn,
                 'gathered_rid_json' : gathered_rid_to_phase_json,
        },
        parameters={},
        dist=distl
    ))

    p_las_fofn_fn =   './1-preads_ovl/las-merge-combine/las_fofn.json'
    hasm_p_ctg_fn    = './3-unzip/1-hasm/p_ctg.fa'
    preads_db_fn     = './1-preads_ovl/build/preads.db'
    preads4falcon_fn = './1-preads_ovl/db2falcon/preads4falcon.fasta'

    wf.addTask(gen_task(
            script=TASK_HASM_SCRIPT,
            inputs={
                'preads_db': preads_db_fn,
                'preads4falcon': preads4falcon_fn,
                'las_fofn': p_las_fofn_fn,
                'rid_to_phase_all': concatenated_rid_to_phase_fn,
                'gathered_rid_json' : gathered_rid_to_phase_json,
            },
            outputs={
                'p_ctg': hasm_p_ctg_fn,
            },
            parameters={},
            dist = Dist(
                job_dict=config['job.step.unzip.hasm'],
                NPROC=1,
            )
    ))

    g2h_all_units_fn = './3-unzip/2-htigs/split/all-units-of-work.json'
    dummy_fn = './3-unzip/2-htigs/split/dummy.sh'
    wf.addTask(gen_task(
            script=TASK_GRAPH_TO_H_TIGS_SPLIT_SCRIPT,
            inputs={
                'falcon_asm_done': falcon_asm_done_fn,
                'preads4falcon': preads4falcon_fn,
                'rid_to_phase_all': concatenated_rid_to_phase_fn,
                'gathered_rid_json': gathered_rid_to_phase_json,
                'p_ctg': hasm_p_ctg_fn,
            },
            outputs={
                'split': g2h_all_units_fn,
                'bash_template': dummy_fn,
            },
            parameters={},
            dist = Dist(
                job_dict=config['job.step.unzip.hasm'],
                NPROC=1,
            )
    ))

    wf.refreshTargets()

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
            dist=Dist(
                NPROC=1,
                job_dict=config['job.step.unzip.hasm'],
                use_tmpdir=False,
            ),
    ))

    combined_ph = "./4-polishing/input/combined_ph.fa"
    combined_ph_fai = "./4-polishing/input/combined_ph.fa.fai"
    combined_eg = "./4-polishing/input/combined_edges.txt"
    readtoctg   = "./4-polishing/input/read2ctg.txt"
    ofastq_fn   = "./4-polishing/input/preamble.fastq"


    wf.addTask(gen_task(
        script=TASK_POLISH_PREAMBLE,
        inputs={
            'P': './3-unzip/all_p_ctg.fa',
            'H': './3-unzip/all_h_ctg.fa',
            'RIDTOPHASE' : './3-unzip/0-phasing/gathered-rid-to-phase/rid_to_phase.all',
            'READNAMELOOKUP' : './3-unzip/readnames/readname_lookup.txt',
            'FQ'             : ifastq_fn,
        },
        outputs={
            'COMBINED' : combined_ph,
            'EDGES' : combined_eg,
            'READTOCTG' : readtoctg,
            "FAI"       : combined_ph_fai,
            "FQO"       : ofastq_fn,
        },
        parameters={},
        dist=Dist(
            NPROC=1,
            job_dict=config['job.step.unzip.hasm'],
            use_tmpdir=False,
        ),
    ))

    wf.refreshTargets()

    collected = dict()

    for ctg in CTGS:
        fn = '4-polishing/temp-unphased/{}/aln.bam'.format(ctg)
        collected['ctg'+ctg] = fn
        wf.addTask(gen_task(
            script=TASK_PLACE_UNPHASED,
            inputs={
                "fa"  : combined_ph,
                "fai" : combined_ph_fai,
                "fq"  : ofastq_fn,
                'readtoctg' : readtoctg,
            },
            outputs={
                "POL": fn,
            },
            parameters={
                'ctg': ctg,
            },
            dist=dist
        ))

    wf.refreshTargets()
    merged_unphased = "4-polishing/merged-unphased/merged_unphased.sorted.bam"

    wf.addTask(gen_task(
        script=TASK_GATHER_UNPHASED,
        inputs=collected,
        outputs={
            "UBAM": merged_unphased,
        },
        parameters={},
        dist=dist
    ))

    wf.refreshTargets()

    SIZES = {}

    with open(combined_ph_fai) as fai:
        for line in fai:
            lc = line.strip().split("\t")
            SIZES[lc[0]] = lc[1]

    PH = {}

    with open(readtoctg) as f:
        for line in f:
            if line.startswith('#'):
                continue
            lc = line.strip().split(" ")
            PH[lc[1]] = 1 + PH.get(lc[1], 0)

    DEST = {}
    for (ctg, count) in PH.iteritems():
        if count < 10:
            LOG.warning("ctg {} is being skipped due to depth of coverage {} < 10 reads total".format(ctg, count))
            DEST[ctg] = 1
        if int(SIZES[ctg]) < 10000:
            LOG.warning("ctg {} is being skipped due to size {} < 10kbp".format(ctg, SIZES[ctg]))
            DEST[ctg] = 1

    for d in DEST:
        del PH[d]

    fns = list()
    LOG.info('len(PH)={}, {!r}'.format(len(PH), dist))
    for ctg in PH:
        fn = "4-polishing/temp-phased/{}/{}.polished.fa".format(ctg, ctg)
        fns.append(fn)
        wf.addTask(gen_task(
            script=TASK_POLISH,
            inputs={
                "FA" : combined_ph,
                "FQ" : ifastq_fn,
                "READTOCTG" : readtoctg,
                "UBAM"      : merged_unphased
            },
            outputs={
                "POL": fn,
            },
            parameters={
                'ctg': ctg,
            },
            dist=dist
        ))

    wf.refreshTargets()

    h_ctg_fn = 'h_ctg.polished.fa'
    p_ctg_fn = 'p_ctg.polished.fa'
    io.touch(h_ctg_fn)
    io.touch(p_ctg_fn)
    for fn in fns:
        if is_haplotig(fn):
            call = 'cat {} >> {}'.format(fn, h_ctg_fn)
        else:
            call = 'cat {} >> {}'.format(fn, p_ctg_fn)
        io.syscall(call)

def is_haplotig(fn):
    """
    >>> is_haplotig('/a/foo_bar.fa')
    True
    >>> is_haplotig('/a/foo.fa')
    False
    """
    return '_' in os.path.basename(fn)

def check(a, b):
    """Simple runtime equality checking.
    """
    if a != b:
        LOG.warning('a != b\n{!r} !=\n{!r}'.format(a, b))
    assert a == b
