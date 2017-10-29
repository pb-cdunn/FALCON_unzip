import argparse
import logging
import sys
from .. import bam_partition_and_merge


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='Partition BAM inputs into BAM files of a small number of ctgs with all the reads for each ctg. These can be easily segregated in parallel later.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--read2ctg-fn', type=str,
        default='./4-quiver/select_reads/read2ctg.msgpack',
        help='Input msgpack from prev step.')
    parser.add_argument(
        '--merged-fn', type=str,
        default='./4-quiver/merge_reads/merged.fofn',
        help='FOFN of merged.bam files. The actual merged files will (probably) be in subdirs of the same directory.')
    parser.add_argument(
        '--max-n-open-files', type=int, default=300,
        help='We write sam files several at-a-time, limited by this.')
    parser.add_argument(
        'input_bam_fofn', type=str,
        help='File of BAM filenames. Paths are relative to dir of FOFN, not CWD.')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(message)s',
    )
    bam_partition_and_merge.run(**vars(args))


if __name__ == '__main__': # pragma: no cover
    main()
