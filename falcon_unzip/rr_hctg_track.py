"""Run in basedir.
"""
from __future__ import print_function, division
from falcon_kit.multiproc import Pool
import falcon_kit.util.io as io
import glob
import msgpack  # for serdes
import os
from heapq import heappush, heappop, heappushpop

Reader = io.CapturedProcessReaderContext

# GLOBALS rid_to_ctg, rid_to_phase


def serialized(rtn):
    return msgpack.dumps(rtn)


def deserialized(rtn_string):
    return msgpack.loads(rtn_string)


def get_rid_to_ctg(fn):
    rid_to_ctg = {}  # local
    with open(fn) as f:
        for row in f:
            row = row.strip().split()
            pid, rid, oid, ctg = row
            rid_to_ctg.setdefault(rid, set())
            rid_to_ctg[rid].add(ctg)
    return rid_to_ctg


def run_tr_stage1(db_fn, fn, min_len, bestn):
    cmd = 'LA4Falcon -m %s %s' % (db_fn, fn)
    reader = Reader(cmd)
    with reader:
        rtn = tr_stage1(reader.readlines, min_len, bestn)
    fn_rtn = '{}.rr_hctg_track.partial.msgpack'.format(fn)
    io.LOG('Ser {!r}'.format(fn_rtn))
    with open(fn_rtn, 'w') as writer:
        writer.write(serialized(rtn))
    return fn_rtn


def tr_stage1(readlines, min_len, bestn):
    """
    for each read in the b-read column inside the LAS files, we
    keep top `bestn` hits with a priority queue through all overlaps
    """
    rtn = {}
    for l in readlines():
        l = l.strip().split()
        q_id, t_id = l[:2]
        overlap_len = -int(l[2])
        idt = float(l[3])
        q_s, q_e, q_l = int(l[5]), int(l[6]), int(l[7])
        t_s, t_e, t_l = int(l[9]), int(l[10]), int(l[11])
        if t_l < min_len:
            continue
        if q_id not in rid_to_ctg:
            continue

        t_phase = rid_to_phase[int(t_id)]
        if t_phase != None:
            ctg_id, block, phase = t_phase
            if block != -1:
                q_phase = rid_to_phase[int(q_id)]
                if q_phase != None:
                    if q_phase[0] == ctg_id and q_phase[1] == block and q_phase[2] != phase:
                        continue

        rtn.setdefault(t_id, [])
        if len(rtn[t_id]) < bestn:
            heappush(rtn[t_id], (overlap_len, q_id))
        else:
            heappushpop(rtn[t_id], (overlap_len, q_id))

    return rtn


def define_global_constants(phased_read_file_fn, read_to_contig_map_fn, rawread_ids_fn):
    global rid_to_ctg, rid_to_phase
    io.LOG('defining constants for track_reads')
    io.logstats()
    rid_to_ctg = get_rid_to_ctg(read_to_contig_map_fn)
    assert rid_to_ctg, 'Empty rid_to_ctg, from {!r}'.format(read_to_contig_map_fn)
    # return here if rid_to_phase not needed

    oid_to_phase = {}
    with open(phased_read_file_fn) as f:
        for row in f:
            row = row.strip().split()
            ctg_id, block, phase = row[1:4]
            oid = row[6]
            block = int(block)
            phase = int(phase)
            oid_to_phase[oid] = (ctg_id, block, phase)
    rid_to_phase = {}
    rid_to_oid = open(rawread_ids_fn).read().split('\n')
    rid_to_phase = [None] * len(rid_to_oid)
    for rid, oid in enumerate(rid_to_oid):
        rid_to_phase[rid] = oid_to_phase.get(oid, None)


def run_track_reads(exe_pool, file_list, min_len, bestn, db_fn):
    io.LOG('running run_track_reads (tr_stage1)')
    io.logstats()
    inputs = []
    io.LOG('num files:{}'.format(len(file_list)))
    for fn in file_list:
        inputs.append((run_tr_stage1, db_fn, fn, min_len, bestn))
    # For each .las input, store the returned dict in a file, named by convention.
    # (See finish_track_reads().)
    fn_rtns = [f for f in exe_pool.imap(io.run_func, inputs)]
    io.LOG('Wrote {} partial.serial files (e.g. {!r}).'.format(len(fn_rtns), fn_rtns[0]))


def finish_track_reads(read_to_contig_map_fn, file_list, bestn, db_fn, rawread_to_contigs_fn):
    io.LOG('running finish_track_reads()')
    rid_to_ctg = get_rid_to_ctg(read_to_contig_map_fn)
    io.LOG('Got rid_to_ctg.')
    # Assume the files already exist, with this naming convention.
    fn_rtns = ['{}.rr_hctg_track.partial.msgpack'.format(fn) for fn in file_list]
    """
    Aggregate hits from each individual LAS and keep the best n hit.
    Note that this does not guarantee that the final results is globally the best n hits espcially
    when the number of `bestn` is too small.  In those case, if there is more hits from single LAS
    file, then we will miss some good  hits.
    """
    bread_to_areads = {}
    for fn_rtn in fn_rtns:
        io.LOG('Deser {!r}'.format(fn_rtn))
        with open(fn_rtn) as reader:
            rtn_string = reader.read()
        res = deserialized(rtn_string)
        del rtn_string
        for k in res:
            bread_to_areads.setdefault(k, [])
            for item in res[k]:
                if len(bread_to_areads[k]) < bestn:
                    heappush(bread_to_areads[k], item)
                else:
                    heappushpop(bread_to_areads[k], item)
        del res

    # rid_to_oid can be helpful for debugging, but otherwise we do not need it.
    #rid_to_oid = open(os.path.join(rawread_dir, 'dump_rawread_ids', 'rawread_ids')).read().split('\n')

    """
    For each b-read, we find the best contig map throgh the b->a->contig map.
    """
    with open(rawread_to_contigs_fn, 'w') as out_f:
        for bread in bread_to_areads:
            ctg_score = {}
            for s, rid in bread_to_areads[bread]:
                if rid not in rid_to_ctg:
                    continue

                ctgs = rid_to_ctg[rid]
                for ctg in ctgs:
                    ctg_score.setdefault(ctg, [0, 0])
                    ctg_score[ctg][0] += -s
                    ctg_score[ctg][1] += 1

            #oid = rid_to_oid[int(bread)]
            ctg_score = list(ctg_score.items())
            ctg_score.sort(key=lambda k: (k[1][0], k[0]))
            rank = 0

            for ctg, score_count in ctg_score:
                if bread in rid_to_ctg and ctg in rid_to_ctg[bread]:
                    in_ctg = 1
                else:
                    in_ctg = 0
                score, count = score_count
                #print(bread, oid, ctg, count, rank, score, in_ctg, file=out_f)
                print(bread, ctg, count, rank, score, in_ctg, file=out_f)
                rank += 1


class TrackReads(object):
    def __init__(self):
        """Scan for 0-rawreads/m*/raw_reads.*.las
        """
        io.LOG('TrackReads.init')
        rawread_dir = os.path.abspath('0-rawreads')
        self.db_fn = os.path.join(rawread_dir, 'raw_reads.db')  # TODO: Another input
        if not os.path.isfile(self.db_fn):
            # It would crash eventually. Actually, it would *not* crash if there are no .las files, but
            # we still want it to crash in that case.
            msg = 'DAZZLER DB {!r} does not exist.'.format(self.db_fn)
            raise Exception(msg)
        # better logic for finding the las files path or move the logic to extern (taking the --fofn option?)
        self.file_list = glob.glob(os.path.join(rawread_dir, 'm*/raw_reads.*.las'))  # TODO: More input
        self.file_list.sort()
        io.LOG('file list: {!r}'.format(self.file_list))

    def try_finish_track_reads(self, read_to_contig_map, bestn, output):
        io.LOG('starting try_finish_track_reads')
        try:
            finish_track_reads(read_to_contig_map, self.file_list, bestn, self.db_fn, output)
            io.LOG('finished finish_track_reads')
        except:
            io.LOG('Exception in finish_track_reads')

    def try_run_track_reads(self, n_core, phased_read_file, read_to_contig_map, rawread_ids, min_len, bestn):
        io.LOG('starting try_run_track_reads')

        n_core = min(n_core, len(self.file_list))

        define_global_constants(phased_read_file, read_to_contig_map, rawread_ids)
        io.LOG('defined global constants')
        io.logstats()

        # We create the Pools *after* we have globals, so we do not need to pass them.
        # (Not much memory saved, but simpler.)
        exe_pool = Pool(n_core)

        try:
            run_track_reads(exe_pool, self.file_list, min_len, bestn, self.db_fn)
            io.LOG('finished track_reads')
        except:
            io.LOG('terminating track_reads workers...')
            exe_pool.terminate()
            raise


def run1(n_core, phased_read_file, read_to_contig_map, rawread_ids, min_len, bestn):
    TrackReads().try_run_track_reads(n_core, phased_read_file, read_to_contig_map, rawread_ids, min_len, bestn)


def run2(read_to_contig_map, bestn, output):
    TrackReads().try_finish_track_reads(read_to_contig_map, bestn, output)
