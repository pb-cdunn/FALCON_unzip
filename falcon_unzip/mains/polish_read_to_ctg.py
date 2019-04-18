"""
This script generates the json that tracks which reads should polish which ctg.
"""

import argparse
import logging
import os
import sys
from collections import defaultdict

def run(lookup_fn, edges_fn, rid_to_phase_fn, out_read_to_ctg_fn):
    nameLookup = {}
    with open(lookup_fn, 'r') as lf:
        for line in lf:
            lc = line.strip().split("\t")
            nameLookup[lc[0]] = lc[1]

    phaseblock2ctg = defaultdict(dict)

    with open(edges_fn, 'r') as mf:
        for line in mf:
            lc = line.strip().split(" ")
            # if lc[5] == -1 it is unphased therefore the read is not assigned to a haplotig
            if(lc[5] == -1):
                continue
            # if lc[0] is a primary, name does not contain an underscore, the reads default to the primary ctg below
            if("_" not in lc[0]):
                continue

            # the rid to phase file does not contain haplotig names so only keep the prefix as part of the tuple key
            ctgname = lc[0].split("_")
            # the key is the primary name, phase block, and phase
            k = (ctgname[0], lc[5], lc[6])
            # multiple assignments of phase blocks are kept as a second level key in the default dict
            phaseblock2ctg[k][lc[0]] = 1


    of = open(out_read_to_ctg_fn, "w")

    of.write("#ccs_id ctg phase_block\n")
    with open(rid_to_phase_fn, 'r') as rctg:
        for line in rctg:
            lc = line.strip().split(" ")
            ctg = lc[1]
            key = (ctg, lc[2], lc[3])
            ls = phaseblock2ctg.get(key, [ctg])
            # for each of the multiple assignments print the read to contig assignment
            for cts in ls:
                of.write("%s %s %s %s %s\n" % (nameLookup[lc[0]], cts, lc[2], lc[0], lc[3]))

class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter): pass

def parse_args(argv):
    description = 'Final stage of assigning reads to contigs'
    epilog = '''
Write read_to_ctgs as
{ccs_id} {ctg} {phase_block}
'''
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=HelpF,
        epilog=epilog,
    )
    parser.add_argument(
        '--lookup-fn', required=True, help='The CCS name ID looup (%%08d <-> CCS)'
    )
    parser.add_argument(
        '--edges-fn',
        required=True,
        help='Combined p & h edges from 2-htg',
    )
    parser.add_argument(
        '--rid-to-phase-fn',
        required=True,
        help='RID to phase file'
    )
    parser.add_argument(
        '--out-read-to-ctg-fn',
        required=True,
        help='Final read assignment for phased reads'
    )

    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
