from .. import io
import os
import sys
from falcon_kit.FastaReader import FastaReader


def run(ctg_id):
    fn = "h_ctg_all.{ctg_id}.fa".format(ctg_id=ctg_id)
    if io.exists_and_not_empty(fn):
        io.syscall("nucmer -mum p_ctg.{ctg_id}.fa h_ctg_all.{ctg_id}.fa -p hp_aln".format(ctg_id=ctg_id))
        io.syscall("show-coords -T -H -l -c hp_aln.delta > hp_aln.coor")
    else:
        print>>sys.stderr, 'No h_ctg_all.{ctg_id}.fa, but that is ok. Continue workflow.\n'
        return 0  # it is ok if there is no h_ctg_all.{ctg_id}.fa, don't want to interupt the workflow

    if not os.path.exists("hp_aln.coor"):
        print>>sys.stderr, 'No "hp_aln.coor". Continuing.\n'
        return 0

    filter_out = set()
    with open("hp_aln.coor") as f:
        for row in f:
            row = row.strip().split()
            q_cov = float(row[10])
            idt = float(row[6])
            if q_cov > 99 and idt > 99.9:
                filter_out.add(row[-1])

    p_ctg_to_phase = {}
    with open("p_ctg_path.%s" % ctg_id) as f:
        for row in f:
            row = row.strip().split()
            b_id, ph_id = (int(row[-2]), int(row[-1]))
            p_ctg_to_phase.setdefault(row[0], {})
            p_ctg_to_phase[row[0]].setdefault((b_id, ph_id), 0)
            p_ctg_to_phase[row[0]][(b_id, ph_id)] += 1

    h_ctg_to_phase = {}
    with open("h_ctg_path.%s" % ctg_id) as f:
        for row in f:
            row = row.strip().split()
            b_id, ph_id = (int(row[-2]), int(row[-1]))
            h_ctg_to_phase.setdefault(row[0], {})
            h_ctg_to_phase[row[0]].setdefault((b_id, ph_id), 0)
            h_ctg_to_phase[row[0]][(b_id, ph_id)] += 1

    with open("h_ctg_ids.%s" % ctg_id, "w") as h_ids:
        with open("h_ctg.%s.fa" % ctg_id, "w") as f:
            h_tig_all = FastaReader("h_ctg_all.%s.fa" % ctg_id)
            for r in h_tig_all:
                p_ctg_phase = p_ctg_to_phase.get(r.name.split("_")[0], {})

                if len(r.sequence) < 500:
                    continue

                if r.name in filter_out:
                    edge_count = sum(h_ctg_to_phase[r.name].values())
                    same_phase_to_p_ctg_count = 0
                    for b_id, ph_id in h_ctg_to_phase[r.name]:
                        if b_id != -1:
                            if (b_id, ph_id) in p_ctg_phase:
                                same_phase_to_p_ctg_count += h_ctg_to_phase[r.name][(b_id, ph_id)]
                    unphased_edge_count = h_ctg_to_phase[r.name] .get((-1, 0), 0)

                    print r.name, edge_count, unphased_edge_count, same_phase_to_p_ctg_count
                    if edge_count - unphased_edge_count - same_phase_to_p_ctg_count < 5:  # there are many non-p_ctg phase segment, do not filter out
                        continue

                print >>f, ">" + r.name
                print >>f, r.sequence
                print >> h_ids, r.name

######
import argparse
import sys


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='Dedup h_tigs for a given ctg_id.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'ctg_id', type=str,
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    run(**vars(args))


if __name__ == "__main__": # pragma: no cover
    main()
