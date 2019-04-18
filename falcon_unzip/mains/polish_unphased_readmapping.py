"""
This script sets up unphased read mapping for polishing.
"""
import argparse
import logging
import os
import sys
import csv
from collections import defaultdict

LOG = logging.getLogger(__name__)

def run(ctg, fai, read_to_ctg, out_ref_names, out_read_names):

    primary2haplotigs = defaultdict(list)

    with open(fai, 'r') as lf:
        for line in lf:
            lc = line.strip().split("\t")
            if '_' in lc[0]:
                ls = lc[0].split("_")
                primary2haplotigs[ls[0]].append(lc[0])

    on = open(out_ref_names, 'w')
    on.write("%s\n" % ctg)

    for seq in primary2haplotigs[ctg]:
        on.write("%s\n" % seq)

    rn = open(out_read_names, 'w')

    with open(read_to_ctg, 'r') as lf:
        # skips header line
        next(lf)
        for line in lf:
            lc = line.strip().split(" ")
            if lc[1] != ctg:
                continue
            if lc[2] != "-1":
                continue

            rn.write("%s\n" % lc[0])


def parse_args(argv):
    description = 'setup for un-phased read tracking'
    parser = argparse.ArgumentParser(
        description=description,
    )
    parser.add_argument(
        '--ctg', required=True, help='primary contig name (p)'
    )
    parser.add_argument(
        '--fai', required=True, help='Fai file of the concated p&h'
    )
    parser.add_argument(
        '--read-to-ctg',
        required=True,
        help='Reads to ctg file that includes the phase blocks',
    )
    parser.add_argument(
        '--out-ref-names',
        required=True,
        help='Sequence names of p&h to map to',
    )
    parser.add_argument(
        '--out-read-names',
        required=True,
        help='Readnames of unphased reads to map to reference',
    )

    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
