from .io import (serialize, deserialize, yield_bam_fn, log, AlignmentFile, mkdirs)
import argparse
import collections
import heapq
import logging
import math
import os
import sys

LOG = logging.getLogger()

def get_bam_header(input_bam_fofn_fn):
    log('Getting BAM header')
    header = None
    for fn in yield_bam_fn(input_bam_fofn_fn):
        with AlignmentFile(fn, 'rb', check_sq = False) as samfile:
            if header is None:
                header = samfile.header
            else:
                header['RG'].extend(samfile.header['RG'])
    try:
        PG = header.pop('PG') #remove PG line as there might be a bug that generates no readable chrs
        #print PG
    except KeyError:
        pass
    log(" Num records in header:", len(header))
    #import pprint
    #log("header:", pprint.pformat(header))
    return header

def partition_ctgs(read2ctg, max_n_groups):
    """
    * Get a count of reads for each ctg, for balanced partitioning.
    * Sort reads, so they might tend to be in the same BAM file.
    * Partition ctgs.
    Return list[set(ctgs)]
    """
    log('Partitioning ctgs into at most {} groups (where num reads is {})'.format(
        max_n_groups, len(read2ctg)))
    log(' Counting reads per ctg')
    ctg2count = collections.defaultdict(int)
    n_reads = 0
    for read, ctg in read2ctg.iteritems():
        ctg2count[ctg] += 1
        n_reads += 1
    log(' {} ctgs, {} reads, => approx. {:.1f} reads per ctg'.format(
        len(ctg2count), len(read2ctg), float(len(read2ctg))/len(ctg2count)))
    assert n_reads == len(read2ctg)

    log('  Sorting reads')
    reads = list(read2ctg.keys())
    reads.sort()

    n_ctgs = len(ctg2count)
    ctgs_per_file = int(math.ceil(n_ctgs / float(max_n_groups)))
    log('  Partitioning {} ctgs; expecting {} per group, given max_n_groups={}'.format(
        n_ctgs, ctgs_per_file, max_n_groups))

    # Sort for highest read-counts at end.
    countedctgs = [(v,k) for (k,v) in ctg2count.items()]
    countedctgs.sort()

    groups = list() # of set(ctgs)
    def yield_groups():
        #nonlocal groups # if needed
        while len(groups) < max_n_groups:
            group = set()
            groups.append(group)
            yield group
        while True:
            # We could use a priority queue for groups, to balance
            # the read-count in each group pretty well.
            # But we want to balance the ctg-count too.
            # So we will add new ctgs backwards and forward
            # through the groups list, repeatedly.
            groups.reverse()
            for group in groups:
                yield group
    next_group = yield_groups().next # to walk back and forth throug list

    while countedctgs:
        # Select biggest remaining ctg (i.e. highest read-count).
        count, ctg = countedctgs.pop()
        group = next_group()
        group.add(ctg)

    # log groups counts
    counts = [len(group) for group in groups]
    log(' group counts: {!r}'.format(counts))

    # log group read-counts
    read_counts = [sum(ctg2count[ctg] for ctg in group) for group in groups]
    log(' group read-counts: {!r}'.format(read_counts))

    return groups

def merge_and_split_alignments(input_bam_fofn_fn, read2ctg, ctg2samfn, samfn2writer):
    """
    For each AlignmentFile in input_bam_fofn_fn,
      for each record in that file,
        find the ctgs from the query_name (based on read2ctg, already selected)
        and append the record to a new samfile based on the ctg name.
    """
    def yield_record_and_ctg():
      """yield (r, ctg)"""
      for i, fn in enumerate(yield_bam_fn(input_bam_fofn_fn)):
        log('BAM input #{:3d}: {!r}'.format(i, fn))
        seen = 0
        used = 0
        with AlignmentFile(fn, 'rb', check_sq=False) as samfile:
          for r in samfile.fetch(until_eof=True):
            seen += 1
            if r.query_name not in read2ctg:
                #print "Missing:", r.query_name
                continue
            used += 1
            ctg = read2ctg[r.query_name]
            yield (r, ctg)
        log(' Saw {} records. Used {}.'.format(seen, used))

    # Actually write. This can take a looooong time for many large BAMs.
    for (r, ctg) in yield_record_and_ctg():
        samfn = ctg2samfn[ctg]
        #log(' Writing to samfn:{!r}'.format(samfn))
        writer = samfn2writer[samfn]
        writer.write(r)

def open_sam_writers(header, sam_fns):
    samfn2writer = dict()
    n = 0
    n_next = 1
    for samfn in sam_fns:
        writer = AlignmentFile(samfn, 'wb', header=header)
        samfn2writer[samfn] = writer

        # log-logging
        n += 1
        if n == n_next:
            log('{} Opening samfile:{!r}'.format(n, samfn))
            n_next = n_next * 2
    return samfn2writer

def close_sam_writers(writers):
    for writer in writers:
        writer.close()

def yield_filenames(ctg_sets, sam_dir):
    log('Choosing SAM filenames in dir {!r}'.format(sam_dir))
    for i, ctgs in enumerate(ctg_sets):
        label = '{:03d}'.format(i)
        sub_dir = os.path.join(sam_dir, label)
        fn = os.path.join(sub_dir, 'merged.bam')
        #log(' SAM fn: {!r}'.format(fn))
        yield fn

def get_ctg2samfn(filenames, ctg_sets):
    assert len(filenames) == len(ctg_sets)
    ctg2samfn = dict()
    for i in range(len(ctg_sets)):
        fn = filenames[i]
        for ctg in ctg_sets[i]:
            ctg2samfn[ctg] = fn
    #import pprint
    #with open('/scratch/cdunn/ctg2samfn.py', 'w') as ofs:
    #    ofs.write(pprint.pformat(ctg2samfn))
    return ctg2samfn

def write_read2ctg_subsets(read2ctg, ctg2samfn):
    """For each samfn, also dump a subset of read2ctg,
    including all the reads which map to any of the ctgs which map to each samfn.
    TODO(CD): Is it possible for a read in read2ctg not to end up the sam-file?
    I think that would mean the BAM inputs do not include all the reads we passed
    to DALIGNER.
    """
    log(" Invert ctg2samfn")
    samfn2ctgs = collections.defaultdict(set)
    for ctg, samfn in ctg2samfn.iteritems():
        samfn2ctgs[samfn].add(ctg)
    log(" Invert read2ctg")
    ctg2reads = collections.defaultdict(set)
    for read, ctg in read2ctg.iteritems():
        ctg2reads[ctg].add(read)
    log(" Write read2ctg subset for each samfn")
    for samfn, ctgs in samfn2ctgs.iteritems():
        read2ctg_subset = dict()
        for ctg in ctgs:
            reads = ctg2reads[ctg]
            for read in reads:
                read2ctg_subset[read] = ctg
        fn = '{}.read2ctg.msgpack'.format(samfn)
        serialize(fn, read2ctg_subset)
        fn = '{}.read2ctg.json'.format(samfn)
        serialize(fn, read2ctg_subset)

def bam_partition_and_merge(input_bam_fofn, read2ctg_fn, merged_fn, max_n_open_files):
    sam_dir = os.path.dirname(merged_fn)
    log('SAM files will go into {!r}'.format(sam_dir))
    mkdirs(sam_dir)
    read2ctg = deserialize(read2ctg_fn)
    ctg_sets = partition_ctgs(read2ctg, max_n_open_files)
    filenames = [fn for fn in yield_filenames(ctg_sets, sam_dir)]
    mkdirs(*[os.path.dirname(fn) for fn in filenames])
    ctg2samfn = get_ctg2samfn(filenames, ctg_sets)
    write_read2ctg_subsets(read2ctg, ctg2samfn)
    header = get_bam_header(input_bam_fofn)
    samfn2writer = open_sam_writers(header, set(ctg2samfn.values()))
    try:
        merge_and_split_alignments(input_bam_fofn, read2ctg, ctg2samfn, samfn2writer)
    finally:
        log('Closing {} SAM writers'.format(len(filenames)))
        close_sam_writers(samfn2writer.values())
    log('Finished partitioning BAM.')
    with open(merged_fn, 'w') as ofs:
        for samfn in samfn2writer.iterkeys():
            ofs.write(samfn + '\n')
    log('Wrote {!r}'.format(merged_fn))

def parse_args(argv):
    parser = argparse.ArgumentParser(description='Partition BAM inputs into BAM files of a small number of ctgs with all the reads for each ctg. These can be easily segregated in parallel later.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--read2ctg-fn', type=str,
            default='./4-quiver/select_reads/read2ctg/read2ctg.msgpack',
            help='Input msgpack from prev step.')
    parser.add_argument('--merged-fn', type=str,
            default='./4-quiver/select_reads/merged/merged.fofn',
            help='FOFN of merged.bam files. The actual merged files will (probably) be in subdirs of the same directory.')
    parser.add_argument('--max-n-open-files', type=int, default=300, help='We write sam files several at-a-time, limited by this.')
    parser.add_argument('input_bam_fofn', type=str, help='File of BAM filenames. Paths are relative to dir of FOFN, not CWD.')
    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)
    bam_partition_and_merge(**vars(args))

if __name__ == '__main__':
    logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s %(message)s',
    )
    main()
