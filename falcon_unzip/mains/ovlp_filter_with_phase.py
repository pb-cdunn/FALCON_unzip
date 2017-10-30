import argparse
import sys
from .. import ovlp_filter_with_phase


def parse_args(argv):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='a simple multi-processes LAS ovelap data filter')
    parser.add_argument(
        '--n_core', type=int, default=4,
        help='number of processes used for generating consensus')
    parser.add_argument(
        '--fofn', type=str, help='file contains the path of all LAS file to be processed in parallel')
    parser.add_argument(
        '--db', type=str, help='read db file path')
    parser.add_argument(
        '--max_diff', type=int, help="max difference of 5' and 3' coverage")
    parser.add_argument(
        '--max_cov', type=int, help="max coverage of 5' or 3' coverage")
    parser.add_argument(
        '--min_cov', type=int, help="min coverage of 5' or 3' coverage")
    parser.add_argument(
        '--min_len', type=int, default=2500, help="min length of the reads")
    parser.add_argument(
        '--bestn', type=int, default=10,
        help="output at least best n overlaps on 5' or 3' ends if possible")
    parser.add_argument(
        '--rid_phase_map', type=str,
        help="the file that encode the relationship of the read id to phase blocks", required=True)
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    ovlp_filter_with_phase.run(args)


if __name__ == '__main__': # pragma: no cover
    main()
