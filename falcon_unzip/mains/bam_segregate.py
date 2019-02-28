from .. import io
from .bam_partition_and_merge import (get_zmw, get_zmw2ctg)
import logging
import os
log = io.log


def segregate_ctgs(merged_fn, zmw2ctg, ctg2samfn, samfn2writer):
    """Walk sequentially through the merged.bam file.
    Write each read to the samfile for its contig.
    (Expensive. Lots of i/o.)
    """
    log('Segregating reads from a merged BAM: {!r}'.format(merged_fn))

    def yield_record_and_ctg():
        """yield (r, ctg)"""
        log('BAM input: {!r}'.format(merged_fn))
        seen = 0
        used = 0
        with io.AlignmentFile(merged_fn, 'rb', check_sq=False) as samfile:
            for r in samfile.fetch(until_eof=True):
                seen += 1
                subread_name = r.query_name # aka "read"
                zmw = get_zmw(subread_name)
                if zmw not in zmw2ctg:
                    # print "Missing:", r.query_name
                    continue
                used += 1
                ctg = zmw2ctg[zmw]
                yield (r, ctg)
        log(' Saw {} records. Used {}.'.format(seen, used))

    # Actually write.
    for (r, ctg) in yield_record_and_ctg():
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


def run(merged_bam_fn, segregated_bam_fns_fn):
    merged_fn = merged_bam_fn
    output_basedir = os.path.normpath(os.path.dirname(segregated_bam_fns_fn))
    if os.path.islink(merged_fn):
        # If it was symlinked, we need to resolve the symlink first, for the implicit dep below.
        # This is kinda hacky. When we make this extra dep explicit, we can drop this hack.
        merged_fn = os.path.realpath(merged_fn)
    # We have (for now) an implicit input next to each merged_fn.
    read2ctg_fn = merged_fn + '.read2ctg.msgpack'  # by convention
    zmw2ctg = get_zmw2ctg(io.deserialize(read2ctg_fn))
    ctg2samfn = get_ctg2samfn(zmw2ctg, output_basedir)
    header = get_single_bam_header(merged_fn)
    bamfns = list(set(ctg2samfn.values()))
    samfn2writer = open_sam_writers(header, bamfns)
    try:
        segregate_ctgs(merged_fn, zmw2ctg, ctg2samfn, samfn2writer)
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
