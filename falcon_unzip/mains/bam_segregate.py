from .. import io
from .bam_partition_and_merge import (get_zmw, get_zmw2ctg)
import logging
import os

LOG = logging.getLogger()


def log(*msgs):
    LOG.info(' '.join(repr(m) for m in msgs))

def yield_record_and_ctg(merged_fn, name2ctg, subread_name2name):
    """yield (r, ctg)

      for each record in the merged BAM file,
        find the ctgs from the query_name
    """
    log('len(name2ctg)=={}'.format(len(name2ctg)))
    log('BAM input: {!r}'.format(merged_fn))
    seen = 0
    used = 0
    with io.AlignmentFile(merged_fn, 'rb', check_sq=False) as samfile:
        for r in samfile.fetch(until_eof=True):
            seen += 1
            subread_name = r.query_name # aka "read"
            name = subread_name2name(subread_name)
            if name not in name2ctg:
                # print "Missing:", r.query_name
                continue
            used += 1
            ctg = name2ctg[name]
            yield (r, ctg)
    log(' Saw {} records. Used {}.'.format(seen, used))

def gen_pairs(merged_fn, read2ctg):
    return yield_record_and_ctg(merged_fn, read2ctg, lambda x: x)

def gen_pairs_extra(merged_fn, read2ctg):
    log('Including all subreads for each mapped zmw.')
    log('len(read2ctg)=={}'.format(len(read2ctg)))
    zmw2ctg = get_zmw2ctg(read2ctg)

    return yield_record_and_ctg(merged_fn, zmw2ctg, get_zmw)

def yield_record_and_ctg_extra(merged_fn, read2ctg):
    """yield (r, ctg)

      for each record in the merged BAM file,
        find the ctgs from the query_name

    Keep all subreads in each mapped zmw.
    """
    name2ctg = get_zmw2ctg(read2ctg)

    log('BAM input: {!r}'.format(merged_fn))
    seen = 0
    used = 0
    with io.AlignmentFile(merged_fn, 'rb', check_sq=False) as samfile:
        for r in samfile.fetch(until_eof=True):
            seen += 1
            subread_name = r.query_name # aka "read"
            name = get_zmw(subread_name)
            if name not in name2ctg:
                # print "Missing:", r.query_name
                continue
            used += 1
            ctg = name2ctg[name]
            yield (r, ctg)
    log(' Saw {} records. Used {}.'.format(seen, used))


def segregate_ctgs(record_and_ctg_pairs, ctg2samfn, samfn2writer):
    """For each (record, ctg) pair,
    Write each read to the samfile for its contig.
    (Expensive. Lots of i/o.)
    """
    # Actually write.
    for (r, ctg) in record_and_ctg_pairs:
        samfn = ctg2samfn[ctg]
        #log(' Writing to samfn:{!r}'.format(samfn))
        writer = samfn2writer[samfn]
        writer.write(r)


def open_sam_writers(header, sam_fns):
    """There will be an open file for each
    sam filename. (Might that be too many?)
    """
    log('Opening {} sam writers'.format(len(sam_fns)))
    io.mkdirs(*[os.path.dirname(fn) for fn in sam_fns])
    samfn2writer = dict()
    n = 0
    n_next = 1
    for samfn in sam_fns:
        # log-logging
        n += 1
        if n == n_next:
            log('{} Opening samfile:{!r}'.format(n, samfn))
            n_next = n_next * 2
        writer = io.AlignmentFile(samfn, 'wb', header=header)
        samfn2writer[samfn] = writer
    return samfn2writer


def close_sam_writers(writers):
    for writer in writers:
        writer.close()


def get_single_bam_header(fn):
    log('Reading BAM header from {!r}'.format(fn))
    with io.AlignmentFile(fn, 'rb', check_sq=False) as samfile:
        return samfile.header


def get_ctg2samfn(zmw2ctg, basedir):
    """Choose the output filename for each BAM-file,
    where each has reads for a single contig.
    """
    ctg2samfn = dict()
    ctgs = set(zmw2ctg.values())
    for ctg in ctgs:
        fn = os.path.join(basedir, '..', 'segregated', ctg, '{}.bam'.format(ctg))
        ctg2samfn[ctg] = fn
    return ctg2samfn


def run(extra_subreads, merged_bam_fn, segregated_bam_fns_fn):
    merged_fn = merged_bam_fn
    output_basedir = os.path.normpath(os.path.dirname(segregated_bam_fns_fn))
    if os.path.islink(merged_fn):
        # If it was symlinked, we need to resolve the symlink first, for the implicit dep below.
        # This is kinda hacky. When we make this extra dep explicit, we can drop this hack.
        merged_fn = os.path.realpath(merged_fn)
    # We have (for now) an implicit input next to each merged_fn.
    read2ctg_fn = merged_fn + '.read2ctg.msgpack'  # by convention
    read2ctg = io.deserialize(read2ctg_fn)
    ctg2samfn = get_ctg2samfn(read2ctg, output_basedir)
    header = get_single_bam_header(merged_fn)
    bamfns = list(set(ctg2samfn.values()))
    samfn2writer = open_sam_writers(header, bamfns)
    try:
        if bool(int(extra_subreads)):
            yield_func = gen_pairs_extra
        else:
            yield_func = gen_pairs
        pairs = yield_func(merged_fn, read2ctg)
        segregate_ctgs(pairs, ctg2samfn, samfn2writer)
    finally:
        close_sam_writers(samfn2writer.values())
    io.serialize(segregated_bam_fns_fn, bamfns)


######
import argparse
import logging
import sys


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='Segregate merged BAM files into single-ctg BAM files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--segregated-bam-fns-fn', type=str,
        help='Output. A JSON list of BAM filenames for segregated reads. The paths must encode each ctg somehow (by convention).',
        # The lines must each have a tuple of fasta and fastq, I think.
    )
    parser.add_argument(
        '--merged-bam-fn', type=str,
        help='Input. A JSON list of (merged BAM, read2ctg).',
    )
    parser.add_argument(
        '--extra-subreads', type=str, default='0',
        help='If set, then include extra subreads (all subreads from any mapped zmws).',
    )
    # parser.add_argument('--max-n-open-files', type=int,
    #        default=300,
    #        help='We write BAM files several at-a-time, hopefully not exceeding this limit.',
    #)
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(message)s',
    )
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
