import argparse
import sys
from .. import graphs_to_h_tigs

def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='layout haplotigs from primary assembly graph and phased aseembly graph')

    parser.add_argument(
        '--fc_asm_path', type=str,
        help='path to the primary Falcon assembly output directory', required=True)
    parser.add_argument(
        '--fc_hasm_path', type=str,
        help='path to the phased Falcon assembly output directory', required=True)
    parser.add_argument(
        '--ctg_id', type=str, help='contig identifier in the bam file', default="all", required=True)
    parser.add_argument(
        '--base_dir', type=str, default="./",
        help='the output base_dir, default to current working directory')
    parser.add_argument(
        '--rid_phase_map', type=str,
        help="path to the file that encode the relationship of the read id to phase blocks", required=True)
    parser.add_argument(
        '--fasta', type=str, help="sequence file of the p-reads", required=True)
    args = parser.parse_args(argv[1:])

    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    graphs_to_h_tigs.run(args)


if __name__ == '__main__': # pragma: no cover
    main()
