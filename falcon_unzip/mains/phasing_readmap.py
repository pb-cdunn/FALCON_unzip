from __future__ import division
import os
import re


def run(out_stream, phased_reads, rawread_ids_fn, pread_ids_fn, pread_to_contigs_fn, the_ctg_id):
    rawread_id_file = rawread_ids_fn
    pread_id_file = pread_ids_fn
    rid_to_oid = open(rawread_id_file).read().split('\n')  # daligner raw read id to the original ids
    pid_to_fid = open(pread_id_file).read().split('\n')  # daligner pread id to the fake ids

    def pid_to_oid(pid):
        fid = pid_to_fid[int(pid)]
        rid = int(fid.split('/')[1]) // 10
        return rid_to_oid[int(rid)]

    rid_to_phase = {}
    with open(phased_reads) as f:
        for row in f:
            row = row.strip().split()
            rid_to_phase[row[6]] = (int(row[2]), int(row[3]))

    arid_to_phase = {}
    map_fn = pread_to_contigs_fn
    with open(map_fn) as f:
        for row in f:
            row = row.strip().split()
            ctg_id = row[1]
            if not ctg_id.startswith(the_ctg_id):
                continue
            if int(row[3]) != 0:  # not the best hit
                continue
            o_id = pid_to_oid(row[0])
            phase = rid_to_phase.get(o_id, (-1, 0))
            arid_to_phase['%09d' % int(row[0])] = phase

    for arid, phase in arid_to_phase.items():
        print >>out_stream, arid, the_ctg_id, phase[0], phase[1]


######
import argparse
import logging
import sys


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='mapping internal daligner read id to phase block and phase',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # we can run this in parallel mode in the furture
    # parser.add_argument('--n_core', type=int, default=4,
    #                    help='number of processes used for generating consensus')
    parser.add_argument(
        '--phased-reads', required=True,
        help='path to read vs. phase map')
    parser.add_argument(
        '--rawread-ids-fn', required=True,
        help='Input. Typically 3-unzip/reads/dump_rawread_ids/rawreads_ids')
    parser.add_argument(
        '--pread-ids-fn', required=True,
        help='Input. Typically 3-unzip/reads/dump_rawread_ids/rawreads_ids')
    parser.add_argument(
        '--pread-to-contigs-fn', required=True,
        help='Input. Typically 3-unzip/reads/pread_to_contigs')
    parser.add_argument(
        '--the-ctg-id', required=True,
        help='contig identifier in the bam file')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig()
    run(sys.stdout, **vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
