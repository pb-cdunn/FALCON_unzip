from .. import io
import os


def get_phased_reads(phased_reads_fn, q_id_map_fn, vmap_fn, p_variant_fn, ctg_id):
    rid_map = io.deserialize(q_id_map_fn)

    read_to_variants = {}
    variant_to_reads = {}
    with open(vmap_fn) as f:
        for l in f:
            if l.startswith('#EOF'): # prove file is complete
                break
            if l.startswith('#'): # skip comments
                continue
            l = l.strip().split()
            variant = "_".join(l[:3])
            read_id = int(l[3])
            read_to_variants.setdefault(read_id, set())
            read_to_variants[read_id].add(variant)
            variant_to_reads.setdefault(variant, set())
            variant_to_reads[variant].add(read_id)
        else:
            raise Exception('No EOF found in {!r}'.format(os.path.abspath(vmap_fn)))

    variant_to_phase = {}
    with open(p_variant_fn) as f:
        for l in f:
            """line format example: V 1 6854 6854_A_A 6854_A_G 6854 22781"""
            l = l.strip().split()
            if l[0] != "V":
                # Skip P lines and comments.
                continue
            pb_id = int(l[1])
            variant_to_phase[l[3]] = (pb_id, 0)
            variant_to_phase[l[4]] = (pb_id, 1)

    with open(phased_reads_fn, "w") as out_f:
        for r in read_to_variants:
            vl = {}
            pl = set()
            for v in list(read_to_variants[r]):
                if v in variant_to_phase:
                    p = variant_to_phase[v]
                    vl[p] = vl.get(p, 0) + 1
                    pl.add(p[0])
            pl = list(pl)
            pl.sort()
            for p in pl:
                if vl.get((p, 0), 0) - vl.get((p, 1), 0) > 1:
                    print >> out_f, r, ctg_id, p, 0, vl.get((p, 0), 0), vl.get((p, 1), 0), rid_map[r]
                elif vl.get((p, 1), 0) - vl.get((p, 0), 0) > 1:
                    print >> out_f, r, ctg_id, p, 1, vl.get((p, 0), 0), vl.get((p, 1), 0), rid_map[r]
    # We were unable to avoid running this particular function in a new directory,
    # so we move the old one out of the way to catch any problems in the workflow graph.
    # TODO: DELETE SOON. (In a couple weeks.)
    old_phased_reads_fn = os.path.join(
        os.path.dirname(phased_reads_fn), '..', 'phased_reads')
    if os.path.exists(old_phased_reads_fn):
        os.rename(old_phased_reads_fn, old_phased_reads_fn + '.bak')


######
import argparse
import sys


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='get phased reads',  # better description?
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--ctg-id', required=True,
    )
    parser.add_argument(
        '--vmap-fn', required=True,
        help='an input'
    )
    parser.add_argument(
        '--p-variant-fn', required=True,
        help='an input'
    )
    parser.add_argument(
        '--q-id-map-fn', required=True,
        help='an input'
    )
    parser.add_argument(
        '--phased-reads-fn', required=True,
        help='an output'
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    get_phased_reads(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
