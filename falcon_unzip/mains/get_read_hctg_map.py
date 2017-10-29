import argparse
import logging
import os
import sys
from .. import get_read_hctg_map

#RAWREAD_DIR = './0-rawreads'
#PREAD_DIR = './1-preads_ovl'
ASM_DIR = './2-asm-falcon'
HASM_DIR = './3-unzip'
QUIVER_DIR = './4-quiver'


def parse_args(argv):
    description = """Generate `read_to_contig_map` that contains the
information from the chain of mapping: (contig id, last col) -> (internal p-read id) -> (internal raw-read id) -> (original read id)
It assumes the inputs are already generated.
"""
    epilog = """You should run this in the base run directory (where "./4-quiver/" lives) if you want to use default arguments."""
    parser = argparse.ArgumentParser(
        description=description, epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--rawread-ids-fn', default=os.path.join(ASM_DIR, 'read_maps/dump_rawread_ids/rawread_ids'),
        help='rawread_ids filename')
    parser.add_argument(
        '--pread-ids-fn', default=os.path.join(ASM_DIR, 'read_maps/dump_pread_ids/pread_ids'),
        help='pread ids filename')
    parser.add_argument(
        '--p-ctg-edges-fn', default=os.path.join(HASM_DIR, 'all_p_ctg_edges'),
        help='primary contig edges filename')
    parser.add_argument(
        '--h-ctg-edges-fn', default=os.path.join(HASM_DIR, 'all_h_ctg_edges'),
        help='haplotype contig edges filename')
    parser.add_argument(
        '--h-ctg-ids-fn', default=os.path.join(HASM_DIR, 'all_h_ctg_ids'),
        help='haplotype contig ids filename')
    parser.add_argument(
        '--read-to-contig-map-fn', default=os.path.join(QUIVER_DIR, 'read_maps/read_to_contig_map'),
        help='OUTPUT: read_to_contig_map filename')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig()
    get_read_hctg_map.run(**vars(args))


if __name__ == '__main__': # pragma: no cover
    main()
