import os
import re


def run(args):
    phased_reads = args.phased_reads
    read_map_dir = args.read_map_dir
    the_ctg_id = args.ctg_id
    base_dir = args.base_dir

    rawread_id_file = os.path.join(read_map_dir, 'dump_rawread_ids', 'rawread_ids')
    pread_id_file = os.path.join(read_map_dir, 'dump_pread_ids', 'pread_ids')
    rid_to_oid = open(rawread_id_file).read().split('\n')  # daligner raw read id to the original ids
    pid_to_fid = open(pread_id_file).read().split('\n')  # daligner pread id to the fake ids

    def pid_to_oid(pid):
        fid = pid_to_fid[int(pid)]
        rid = int(fid.split('/')[1]) / 10
        return rid_to_oid[int(rid)]

    rid_to_oid = open(rawread_id_file).read().split('\n')  # daligner raw read id to the original ids
    pid_to_fid = open(pread_id_file).read().split('\n')  # daligner pread id to the fake ids

    rid_to_phase = {}
    with open(phased_reads) as f:
        for row in f:
            row = row.strip().split()
            rid_to_phase[row[6]] = (int(row[2]), int(row[3]))

    arid_to_phase = {}
    map_fn = os.path.join(read_map_dir, 'pread_to_contigs')
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

    with open(os.path.join(base_dir, 'rid_to_phase.%s' % the_ctg_id), 'w') as f:
        for arid, phase in arid_to_phase.items():
            print >>f, arid, the_ctg_id, phase[0], phase[1]
