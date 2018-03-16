#! /usr/bin/env python2.7

import os
import sys
import argparse
import execute
from falcon_kit.FastaReader import FastaReader

def load_rid_to_phase(rid_phase_map, ctg_id):
    arid2phase = {}
    with open(rid_phase_map) as f:
        for row in f:
            row = row.strip().split()
            pread_id, phase_ctg_id, phase_blk_id, phase_id = row[0:4]
            if phase_ctg_id == ctg_id:
                arid2phase[pread_id] = (phase_ctg_id, phase_blk_id, phase_id)
    return arid2phase

def extract_preads(preads_paths, header_set, out_path):
    """
    This is a whitelist filter.
    Takes a set of headers (header_set) which should
    be extracted from the preads FASTA file (preads_path),
    and written to the out_path.
    """
    num_seqs = 0
    fp_preads = FastaReader(preads_paths)
    with open(out_path, 'w') as fp_out:
        for record in fp_preads:
            header = record.name.split()[0]
            seq = record.sequence.upper()
            if header in header_set:
                num_seqs += 1
                fp_out.write('>%s\n%s\n' % (header, seq))

def run(ctg_id, rid_phase_map, preads, out):
    rid2phase = load_rid_to_phase(rid_phase_map, ctg_id)
    pread_header_set = set(rid2phase.keys())
    extract_preads(preads, pread_header_set, out)

def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='This script prepares all info needed for prototyping the Unzip workflow.',  # better description?
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--ctg-id', required=True,
        help='The contig ID to process.'
    )
    parser.add_argument(
        '--rid-phase-map', required=True, type=str,
        help='The file that encodes the relationship of the read id to phase blocks.'
    )
    parser.add_argument(
        '--preads', required=True, type=str,
        help='Path to the preads4falcon.fasta file.'
    )
    parser.add_argument(
        '--out', required=True,
        default='./',
        help='Output file to write the selected preads.',
    )

    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    # This file contains all phased reads: --rid_phase_map 3-unzip/1-hasm/concatenated-rid-to-phase/rid_to_phase.all
    # This file contains all the phasing blocks (lines beginning with 'P'): --p-variant-fn 3-unzip/0-phasing/000000F/get_phased_blocks/phased_variants

    args = parse_args(argv)
    run(**vars(args))

if __name__ == '__main__':  # pragma: no cover
    main()
