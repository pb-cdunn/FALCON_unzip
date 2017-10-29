import argparse
import sys
from .. import phased_ovlp_to_graph


def parse_args(argv):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='a string graph assembler that is desinged for handling diploid genomes')
    parser.add_argument(
        'overlap_file', help='a file that contains the overlap information.')
    parser.add_argument(
        '--min_len', type=int, default=4000,
        help='minimum length of the reads to be considered for assembling')
    parser.add_argument(
        '--min_idt', type=float, default=96,
        help='minimum alignment identity of the reads to be considered for assembling')
    parser.add_argument(
        '--lfc', action="store_true", default=False,
        help='use local flow constraint method rather than best overlap method to resolve knots in string graph')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    phased_ovlp_to_graph.run(args)


if __name__ == '__main__': # pragma: no cover
    main()
