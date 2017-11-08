import os

def generate_association_table(vmap_fn, atable_fn, ctg_id):
    vmap = {}
    v_positions = []

    with open(vmap_fn) as f:
        for l in f:
            if l.startswith('#EOF'): # prove file is complete
                break
            if l.startswith('#'): # skip comments
                continue
            l = l.strip().split()
            pos = int(l[0])
            ref_b = l[1]
            v_b = l[2]
            q_id = int(l[3])
            if (pos, ref_b) not in vmap:
                v_positions.append((pos, ref_b))
            vmap.setdefault((pos, ref_b), {})
            vmap[(pos, ref_b)].setdefault(v_b, [])
            vmap[(pos, ref_b)][v_b].append(q_id)
        else:
            raise Exception('No EOF found in {!r}'.format(os.path.abspath(vmap_fn)))

    #xary = []
    #yary = []
    with open(atable_fn, "w") as out_f:
        for i1 in xrange(len(v_positions)):
            link_count = 0
            for i2 in xrange(i1 + 1, len(v_positions)):
                pos1, rb1 = v_positions[i1]
                pos2, rb2 = v_positions[i2]
                if pos2 - pos1 > (1 << 16):
                    continue
                ct = {}
                p1table = []
                p2table = []
                s1 = 0
                list1 = vmap[(pos1, rb1)].items()
                assert len(list1) == 2, 'len={}'.format(len(list1))
                for b1, qids1 in list1:
                    p1table.append((b1, len(qids1)))
                    s1 += len(qids1)

                s2 = 0
                list2 = vmap[(pos2, rb2)].items()
                assert len(list2) == 2, 'len={}'.format(len(list2))
                for b2, qids2 in list2:
                    p2table.append((b2, len(qids2)))
                    s2 += len(qids2)

                total_s = 0
                for b1, qids1 in list1:
                    for b2, qids2 in list2:
                        s = len(set(qids1) & set(qids2))
                        ct[(b1, b2)] = s
                        total_s += s
                if total_s < 6:
                    continue

                b11 = p1table[0][0]
                b12 = p1table[1][0]
                b21 = p2table[0][0]
                b22 = p2table[1][0]
                print >> out_f, pos1, b11, b12, pos2, b21, b22, ct[(
                    b11, b21)], ct[(b11, b22)], ct[(b12, b21)], ct[(b12, b22)]

                # xary.append(pos1)
                # yary.append(pos2)
                link_count += 1
                if link_count > 500:
                    break


######
import argparse
import sys


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='generate association table',  # better description?
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--ctg-id', required=True,
    )
    parser.add_argument(
        '--vmap-fn', required=True,
        help='an input'
    )
    parser.add_argument(
        '--atable-fn', required=True,
        help='an output'
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    generate_association_table(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
