import argparse
import logging
import sys
from .. import get_read2ctg


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='Map ctg->BAM filename.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--output', type=str,
        default='./4-quiver/select_reads/read2ctg.msgpack',
        help='Serialized map of ctg to list of BAM which contains it.')
    parser.add_argument(
        '--rawread-to-contigs', type=str,
        default='./2-asm-falcon/read_maps/dump_rawread_ids/rawread_to_contigs',
        help='rawread_to_contigs file (from where?)')
    parser.add_argument(
        '--rawread-ids', type=str,
        default='./2-asm-falcon/read_maps/dump_rawread_ids/rawread_ids', help='rawread_ids file (from where?)')
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
    get_read2ctg.write_read2ctg(**vars(args))


if __name__ == "__main__":
    main()
