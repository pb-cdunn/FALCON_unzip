import argparse
import logging
import os
import sys
from pypeflow.simple_pwatcher_bridge import (PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase,
                                             makePypeLocalFile, fn, PypeTask)
from .. import phasing
from ..tasks import unzip as tasks_unzip

LOG = logging.getLogger(__name__)


def get_phasing_tasks(bam_fn, fasta_fn, ctg_id, base_dir, samtools):
    LOG.debug('IN PHASING')

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
    )(tasks_unzip.task_make_het_call)

    yield make_het_call_task

    atable_file = makePypeLocalFile(os.path.join(base_dir, ctg_id, 'g_atable', "atable"))
    parameters = {}
    parameters["ctg_id"] = ctg_id
    parameters["base_dir"] = base_dir
    generate_association_table_task = PypeTask(inputs={"vmap_file": vmap_file},
                                               outputs={"atable_file": atable_file},
                                               parameters=parameters,
                                               )(tasks_unzip.task_generate_association_table)

    yield generate_association_table_task

    phased_variant_file = makePypeLocalFile(os.path.join(base_dir, ctg_id, 'get_phased_blocks', "phased_variants"))
    get_phased_blocks_task = PypeTask(inputs={"vmap_file": vmap_file, "atable_file": atable_file},
                                      outputs={"phased_variant_file": phased_variant_file},
                                      )(tasks_unzip.task_get_phased_blocks)
    yield get_phased_blocks_task

    phased_read_file = makePypeLocalFile(os.path.join(base_dir, ctg_id, "phased_reads"))
    get_phased_reads_task = PypeTask(inputs={"vmap_file": vmap_file,
                                             "q_id_map_file": q_id_map_file,
                                             "phased_variant_file": phased_variant_file},
                                     outputs={"phased_read_file": phased_read_file},
                                     parameters={"ctg_id": ctg_id},
                                     )(tasks_unzip.task_get_phased_reads)
    yield get_phased_reads_task


def run_phasing(args):
    kwds = dict(
        bam_fn = args.bam,
        fasta_fn = args.fasta,
        ctg_id = args.ctg_id,
        base_dir = args.base_dir,
        samtools = args.samtools,
    )
    tasks = list(get_phasing_tasks(**kwds))
    wf = PypeProcWatcherWorkflow(
        max_jobs=1,
    )
    wf.addTasks(tasks)
    wf.refreshTargets()
    # with open("fc_phasing_wf.dot", "w") as f:
    #    print >>f, wf.graphvizDot


def parse_args(argv):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='phasing variants and reads from a bam file')
    # we can run this in parallel mode in the furture
    # parser.add_argument('--n_core', type=int, default=4,
    #                    help='number of processes used for generating consensus')
    parser.add_argument(
        '--bam', type=str, help='path to sorted bam file', required=True)
    parser.add_argument(
        '--fasta', type=str, help='path to the fasta file of contain the contig', required=True)
    parser.add_argument(
        '--ctg_id', type=str, help='contig identifier in the bam file', required=True)
    parser.add_argument(
        '--base_dir', type=str, default="./",
        help='the output base_dir, default to current working directory')
    parser.add_argument(
        '--samtools', type=str, default="samtools", help='path to samtools')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    run_phasing(args)


if __name__ == '__main__': # pragma: no cover
    main()
