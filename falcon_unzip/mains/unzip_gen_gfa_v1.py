import argparse
import sys
from .. import unzip_gen_gfa_v1


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Generates GFA output (on stdout) from FALCON's assembly.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--preads-fasta', type=str, default='preads4falcon.fasta',
        help='path to the preads4falcon.fasta file')
    parser.add_argument(
        '--p-ctg-fasta', type=str, default='all_p_ctg.fa', help='path to the primary contigs file')
    parser.add_argument(
        '--h-ctg-fasta', type=str, default='all_h_ctg.fa', help='path to the haplotig file')
    parser.add_argument(
        '--unzip-root', type=str, default='./', help='path to the 3-unzip directory')
    parser.add_argument(
        '--add-string-graph', action='store_true',
        help="in addition to tiling paths, output other edges and nodes from the final string graph")
    parser.add_argument(
        '--write-reads', '-r', action='store_true', help="output read sequences in S lines")
    parser.add_argument(
        '--write-contigs', '-c', action='store_true', help="output contig sequences as S lines")
    parser.add_argument(
        '--min-p-len', type=int, default=0,
        help='primary contig paths with length smaller than this will not be reported')
    parser.add_argument(
        '--min-h-len', type=int, default=0,
        help='haplotig paths with length smaller than this will not be reported')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    unzip_gen_gfa_v1.gfa_from_unzip(sys.stdout, **vars(args))


if __name__ == '__main__': # pragma: no cover
    main()
