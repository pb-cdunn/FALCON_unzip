from __future__ import division
import os


def make_dirs(d):
    if not os.path.isdir(d):
        os.makedirs(d)


def run(
        base_dir,
        rawread_ids_fn, pread_ids_fn,
        h_ctg_ids_fn,
        h_ctg_edges_fn,
        p_ctg_edges_fn,
        output_fn,
    ):
    #make_dirs(os.path.dirname(output_fn))  # Workflow does this too.

    pread_did_to_rid = open(pread_ids_fn).read().split('\n')
    rid_to_oid = open(rawread_ids_fn).read().split('\n')

    h_ctg_ids = set()
    with open(h_ctg_ids_fn) as f:
        for row in f:
            row = row.strip()
            h_ctg_ids.add(row)

    pread_to_contigs = {}

    for fnname in (p_ctg_edges_fn, h_ctg_edges_fn):
        with open(fnname) as f:
            for row in f:
                row = row.strip().split()
                ctg = row[0]
                if len(ctg.split('_')) > 1 and ctg not in h_ctg_ids:
                    continue
                n1 = row[1]
                n2 = row[2]
                pid1 = int(n1.split(':')[0])
                pid2 = int(n2.split(':')[0])
                rid1 = pread_did_to_rid[pid1].split('/')[1]
                rid2 = pread_did_to_rid[pid2].split('/')[1]
                rid1 = int(int(rid1) // 10)
                rid2 = int(int(rid2) // 10)
                oid1 = rid_to_oid[rid1]
                oid2 = rid_to_oid[rid2]
                k1 = (pid1, rid1, oid1)
                pread_to_contigs.setdefault(k1, set())
                pread_to_contigs[k1].add(ctg)
                k2 = (pid2, rid2, oid2)
                pread_to_contigs.setdefault(k2, set())
                pread_to_contigs[k2].add(ctg)

    assert pread_to_contigs, 'Empty p/h_ctg_edges: {!r} {!r}'.format(
        p_ctg_edges_fn, h_ctg_edges_fn)
    with open(output_fn, 'w') as f:
        for k in pread_to_contigs:
            pid, rid, oid = k
            for ctg in list(pread_to_contigs[k]):
                print >>f, '%09d %09d %s %s' % (pid, rid, oid, ctg)


######
import argparse
import logging
import os
import sys


def parse_args(argv):
    description = """Generate `read_to_contig_map` that contains the
information from the chain of mapping: (contig id, last col) -> (internal p-read id) -> (internal raw-read id) -> (original read id)
It assumes the inputs are already generated.
"""
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description, epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--base-dir', default='.',
        help='Substituted as {base_dir} into default inputs.')
    parser.add_argument(
        '--rawread-ids-fn', default='{base_dir}/3-unzip/reads/dump_rawread_ids/rawread_ids',
        help='rawread_ids filename')
    parser.add_argument(
        '--pread-ids-fn', default='{base_dir}/3-unzip/reads/dump_pread_ids/pread_ids',
        help='pread ids filename')
    parser.add_argument(
        '--p-ctg-edges-fn', default='{base_dir}/3-unzip/all_p_ctg_edges',
        help='primary contig edges filename')
    parser.add_argument(
        '--h-ctg-edges-fn', default='{base_dir}/3-unzip/all_h_ctg_edges',
        help='haplotype contig edges filename')
    parser.add_argument(
        '--h-ctg-ids-fn', default='{base_dir}/3-unzip/all_h_ctg_ids',
        help='haplotype contig ids filename')
    parser.add_argument(
        '--output-fn', default='{base_dir}/4-polish/read_maps/read_to_contig_map',
        help='output read_to_contig_map filename')
    args = parser.parse_args(argv[1:])
    args.rawread_ids_fn = args.rawread_ids_fn.format(**vars(args))
    args.pread_ids_fn = args.pread_ids_fn.format(**vars(args))
    args.h_ctg_ids_fn = args.h_ctg_ids_fn.format(**vars(args))
    args.h_ctg_edges_fn = args.h_ctg_edges_fn.format(**vars(args))
    args.p_ctg_edges_fn = args.p_ctg_edges_fn.format(**vars(args))
    args.output_fn = args.output_fn.format(**vars(args))
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig()
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
