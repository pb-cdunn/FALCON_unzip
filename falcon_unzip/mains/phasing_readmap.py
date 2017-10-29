import argparse
import logging
import sys
from .. import phasing_readmap


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='mapping internal daligner read id to phase block and phase',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # we can run this in parallel mode in the furture
    # parser.add_argument('--n_core', type=int, default=4,
    #                    help='number of processes used for generating consensus')
    parser.add_argument(
        '--phased_reads', type=str, help='path to read vs. phase map', required=True)
    parser.add_argument(
        '--read_map_dir', type=str, help='path to the read map directory', required=True)
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
    phasing_readmap.run(args)


if __name__ == '__main__': # pragma: no cover
    main()
