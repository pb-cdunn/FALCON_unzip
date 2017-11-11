"""Deprecated.
These tasks are now run as part of the full workflow. But in theory, they can
be run separately here.
"""
import argparse
import logging
import os
import sys
from pypeflow.simple_pwatcher_bridge import (PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase,
                                             makePypeLocalFile, fn, PypeTask)
from ..tasks import unzip as tasks_unzip


def run_phasing(args):
    phased_reads_file = makePypeLocalFile(os.path.join(
        args.base_dir, args.ctg_id, 'phased_reads'))
    bam_file = makePypeLocalFile(args.bam)
    fasta_file = makePypeLocalFile(args.fasta)
    kwds = dict(
        phased_reads_file=phased_reads_file,
        bam_file=bam_file,
        fasta_file=fasta_file,
        ctg_id=args.ctg_id,
        base_dir=args.base_dir,
    )
    tasks = list(tasks_unzip.yield_phasing_tasks(**kwds))
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
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    run_phasing(args)


if __name__ == '__main__':  # pragma: no cover
    main()
