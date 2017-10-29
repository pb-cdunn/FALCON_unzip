import argparse
import logging
import sys
from .. import bam_segregate


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='Segregate merged BAM files into single-ctg BAM files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--output-fn', type=str,
        help='A FOFN of BAM for segregated reads. The paths must encode each ctg somehow (by convention).',
    )
    parser.add_argument(
        '--merged-fn', type=str,
        help='A merged BAM file.',
    )
    # parser.add_argument('--max-n-open-files', type=int,
    #        default=300,
    #        help='We write BAM files several at-a-time, hopefully not exceeding this limit.',
    #)
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(message)s',
    )
    bam_segregate.run(**vars(args))


if __name__ == '__main__': # pragma: no cover
    main()
