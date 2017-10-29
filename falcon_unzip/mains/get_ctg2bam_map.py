"""WORK IN PROGRESS

And not currently used. Since ctgs tend to be spread over many BAM files,
this might not be useful information.
"""
from ..io import (serialize, deserialize, yield_bam_fn, log, AlignmentFile)
import argparse
import logging
import sys


def get_ctg2bam(input_bam_fofn_fn, read2ctg):
    return []


def write_ctg2bam(output, input_bam_fofn, int_file):
    read2ctg = deserialize(int_file)
    ctg2bam = get_ctg2bam(input_bam_fofn_fn=input_bam_fofn, read2ctg=read2ctg)
    serialize(output, ctg2bam)
    serialize(output + '.json', ctg2bam)


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='Map ctg->BAM filename.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--int-file', type=str, default='/scratch/cdunn/read2ctg0.msgpack',
        help='Input msgpack from prev step.')
    parser.add_argument(
        '--output', type=str, default='/scratch/cdunn/ctg2bam.msgpack',
        help='Serialized map of ctg to list of BAM which contains it.')
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
    write_ctg2bam(**vars(args))


if __name__ == "__main__": # pragma: no cover
    main()
