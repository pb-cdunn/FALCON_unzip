import os

def get_score(c_score, pos1, pos2, s1, s2):
    if pos1 > pos2:
        pos1, pos2 = pos2, pos1
        s1, s2 = s2, s1
    b11, b12 = s1
    b21, b22 = s2
    return c_score[(pos1, pos2)][(b11 + b21, b12 + b22)]


def get_phased_blocks(vmap_fn, atable_fn, p_variant_fn):
    left_connect = {}
    right_connect = {}

    c_score = {}
    states = {}
    positions = set()

    ref_base = {}
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
            ref_base[pos] = ref_b
        else:
            raise Exception('No EOF found in {!r}'.format(os.path.abspath(vmap_fn)))

    with open(atable_fn) as f:
        for l in f:
            l = l.strip().split()
            pos1, b11, b12, pos2, b21, b22, s11, s12, s21, s22 = l
            s11, s12, s21, s22 = int(s11), int(s12), int(s21), int(s22)
            if abs(s11 + s22 - s12 - s21) < 6:
                continue
            pos1 = int(pos1)
            pos2 = int(pos2)
            positions.add(pos1)
            positions.add(pos2)
            right_connect.setdefault(pos1, [])
            right_connect[pos1].append(pos2)
            left_connect.setdefault(pos2, [])
            left_connect[pos2].append(pos1)
            c_score[(pos1, pos2)] = {(b11 + b21, b12 + b22): s11 + s22, (b12 + b22, b11 + b21): s11 + s22,
                                     (b12 + b21, b11 + b22): s12 + s21, (b11 + b22, b12 + b21): s12 + s21}

            if pos1 not in states:
                st1 = (b11, b12)
                st2 = (b12, b11)
                score1 = 0
                score2 = 0
                for pp in left_connect.get(pos1, []):
                    if pp in states:
                        st0 = states[pp]
                    else:
                        continue
                    score1 += get_score(c_score, pp, pos1, st0, st1)
                    score2 += get_score(c_score, pp, pos1, st0, st2)

                for pp in right_connect.get(pos1, []):
                    if pp in states:
                        st0 = states[pp]
                    else:
                        continue
                    score1 += get_score(c_score, pos1, pp, st1, st0)
                    score2 += get_score(c_score, pos1, pp, st2, st0)

                if score1 >= score2:
                    states[pos1] = st1
                else:
                    states[pos1] = st2

            if pos2 not in states:
                st1 = (b21, b22)
                st2 = (b22, b21)
                score1 = 0
                score2 = 0
                for pp in left_connect.get(pos2, []):
                    if pp in states:
                        st0 = states[pp]
                    else:
                        continue
                    score1 += get_score(c_score, pp, pos2, st0, st1)
                    score2 += get_score(c_score, pp, pos2, st0, st2)

                for pp in right_connect.get(pos2, []):
                    if pp in states:
                        st0 = states[pp]
                    else:
                        continue
                    score1 += get_score(c_score, pos2, pp, st1, st0)
                    score2 += get_score(c_score, pos2, pp, st2, st0)

                if score1 >= score2:
                    states[pos2] = st1
                else:
                    states[pos2] = st2

    positions = list(positions)
    positions.sort()

    iter_count = 0
    while 1:
        iter_count += 1
        if iter_count > 10:
            break
        update_count = 0
        for p in positions:
            b1, b2 = states[p]
            st1 = (b1, b2)
            st2 = (b2, b1)

            score1 = 0
            score2 = 0
            for pp in left_connect.get(p, []):
                st0 = states[pp]
                score1 += get_score(c_score, pp, p, st0, st1)
                score2 += get_score(c_score, pp, p, st0, st2)

            # for pp in right_connect.get(p,[]):
            #    st0 = states[pp]
            #    score1 += get_score( c_score, p, pp, st1 ,st0)
            #    score2 += get_score( c_score, p, pp, st2, st0)

            if score1 >= score2:
                states[p] = st1
            else:
                states[p] = st2
                update_count += 1
        if update_count == 0:
            break

    right_extent = {}
    right_score = {}
    left_extent = {}
    left_score = {}

    for p in positions:

        left_extent[p] = p
        left_score[p] = 0
        if p in left_connect:
            left = p
            st0 = states[p]
            st0_ = st0[1], st0[0]
            for pp in left_connect[p]:
                st1 = states[pp]
                s = get_score(c_score, pp, p, st1, st0)
                s_ = get_score(c_score, pp, p, st1, st0_)
                left_score[p] += s - s_
                if s - s_ > 0 and pp < left:
                    left = pp
            left_extent[p] = left

        right_extent[p] = p
        right_score[p] = 0
        if p in right_connect:
            right = p
            st0 = states[p]
            st0_ = st0[1], st0[0]
            for pp in right_connect[p]:
                st1 = states[pp]
                s = get_score(c_score, p, pp, st0, st1)
                s_ = get_score(c_score, p, pp, st0_, st1)
                right_score[p] += s - s_
                if s - s_ > 0 and pp > right:
                    right = pp
            right_extent[p] = right

    phase_block_id = 1
    phase_blocks = {}
    pb = []

    max_right_ext = 0
    for p in positions:
        if right_score[p] < 10 or left_score[p] < 10:
            continue
        b1, b2 = states[p]
        if max_right_ext < left_extent[p]:
            if len(pb) > 3:
                phase_blocks[phase_block_id] = pb
                phase_block_id += 1
            pb = []
        pb.append((p, b1, b2))
        if right_extent[p] > max_right_ext:
            max_right_ext = right_extent[p]
    if len(pb) > 3:
        phase_blocks[phase_block_id] = pb
    else:
        phase_block_id -= 1

    with open(p_variant_fn, "w") as out_f:
        for pid in xrange(1, phase_block_id + 1):
            if len(phase_blocks[pid]) == 0:
                continue
            min_ = min([x[0] for x in phase_blocks[pid]])
            max_ = max([x[0] for x in phase_blocks[pid]])

            print >>out_f, "P", pid, min_, max_, max_ - \
                min_, len(phase_blocks[pid]), 1.0 * (max_ - min_) / len(phase_blocks[pid])
            for p, b1, b2 in phase_blocks[pid]:
                rb = ref_base[p]
                print >>out_f, "V", pid, p, "%d_%s_%s" % (p, rb, b1), "%d_%s_%s" % (
                    p, rb, b2), left_extent[p], right_extent[p], left_score[p], right_score[p]


######
import argparse
import sys


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='generate association table',  # better description?
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--vmap-fn', required=True,
        help='an input'
    )
    parser.add_argument(
        '--atable-fn', required=True,
        help='an input'
    )
    parser.add_argument(
        '--p-variant-fn', required=True,
        help='an output'
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    get_phased_blocks(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
