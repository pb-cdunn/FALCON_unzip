import argparse
import logging
import os
import sys

def make_dirs(d):
    if not os.path.isdir(d):
        os.makedirs(d)

def generate_read_to_hctg_map(
        rawread_ids_fn, pread_ids_fn,
        h_ctg_ids_fn,
        h_ctg_edges_fn,
        p_ctg_edges_fn,
        read_to_contig_map_fn):
    make_dirs(os.path.dirname(read_to_contig_map_fn)) # Workflow does this too.

    pread_did_to_rid = open(pread_ids_fn).read().split('\n')
    rid_to_oid = open(rawread_ids_fn).read().split('\n')

    h_ctg_ids = set()
    with open(h_ctg_ids_fn) as f:
        for row in f:
            row = row.strip()
            h_ctg_ids.add( row )

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
                rid1 = int(int(rid1)/10)
                rid2 = int(int(rid2)/10)
                oid1 = rid_to_oid[rid1]
                oid2 = rid_to_oid[rid2]
                k1 = (pid1, rid1, oid1)
                pread_to_contigs.setdefault( k1, set() )
                pread_to_contigs[ k1 ].add( ctg )
                k2 = (pid2, rid2, oid2)
                pread_to_contigs.setdefault( k2, set() )
                pread_to_contigs[ k2 ].add( ctg )

    assert pread_to_contigs, 'Empty p/h_ctg_edges: {!r} {!r}'.format(
            p_ctg_edges_fn, h_ctg_edges_fn)
    with open(read_to_contig_map_fn, 'w') as f:
        for k in pread_to_contigs:
            pid, rid, oid = k
            for ctg in list(pread_to_contigs[ k ]):
                print >>f, '%09d %09d %s %s' % (pid, rid, oid, ctg)


#RAWREAD_DIR = './0-rawreads'
#PREAD_DIR = './1-preads_ovl'
ASM_DIR = './2-asm-falcon'
HASM_DIR = './3-unzip'
QUIVER_DIR = './4-quiver'

def parse_args(argv):
    description="""Generate `read_to_contig_map` that contains the
information from the chain of mapping: (contig id, last col) -> (internal p-read id) -> (internal raw-read id) -> (original read id)
It assumes the inputs are already generated.
"""
    epilog = """You should run this in the base run directory (where "./4-quiver/" lives) if you want to use default arguments."""
    parser = argparse.ArgumentParser(description=description, epilog=epilog,
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--rawread-ids-fn', default=os.path.join(ASM_DIR, 'read_maps/dump_rawread_ids/rawread_ids'),
            help='rawread_ids filename')
    parser.add_argument('--pread-ids-fn', default=os.path.join(ASM_DIR, 'read_maps/dump_pread_ids/pread_ids'),
            help='pread ids filename')
    parser.add_argument('--p-ctg-edges-fn', default=os.path.join(HASM_DIR, 'all_p_ctg_edges'),
            help='primary contig edges filename')
    parser.add_argument('--h-ctg-edges-fn', default=os.path.join(HASM_DIR, 'all_h_ctg_edges'),
            help='haplotype contig edges filename')
    parser.add_argument('--h-ctg-ids-fn', default=os.path.join(HASM_DIR, 'all_h_ctg_ids'),
            help='haplotype contig ids filename')
    parser.add_argument('--read-to-contig-map-fn', default=os.path.join(QUIVER_DIR, 'read_maps/read_to_contig_map'),
            help='OUTPUT: read_to_contig_map filename')
    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    logging.basicConfig()

    args = parse_args(argv)

    generate_read_to_hctg_map(**vars(args))
