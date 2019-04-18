"""
This script parses an FAI file and generates a list of primary and assoicate contigs
"""
import argparse
import logging
import os
import sys
import csv

def run(fai, ctg):
    """Reads FAI and prints matching ctgs."""
    with open(fai, 'r') as lf:
        for line in lf:
            lc = line.strip().split("\t")
            name = lc[0].split("_")
            if(name[0] == ctg):
                print(lc[0])


def parse_args(argv):
    description = 'Generates a batch (primary and haplotigs)'
    parser = argparse.ArgumentParser(
        description=description,
    )
    parser.add_argument(
        '--fai', required=True, help='fai for combined p and h contig fastas'
    )
    parser.add_argument(
        '--ctg',
        required=True,
        help='The base name of the contig. E.G. 000001F',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
