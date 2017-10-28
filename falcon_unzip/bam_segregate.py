from . import io
import argparse
import logging
import os
import re
import sys
#LOG = logging.getLogger()
log = io.log

re_merged_fn = re.compile(r'([^/.]+).bam$')


def get_ctg_from_fn(fn):
    """
    >>> get_ctg_from_fn('/a/b/c.bam')
    'c'
    >>> get_ctg_from_fn('/a/b/c12.bam')
    'c12'
    """
    log(' get_ctg_from_fn({!r})'.format(fn))
    mo = re_merged_fn.search(fn)
    if not mo:
        raise Exception('No ctg found in BAM filename {!r} based on regex {!r}'.format(
            fn, re_merged_fn.pattern))
    return mo.group(1)


def segregate_ctgs(merged_fn, read2ctg, ctg2samfn, samfn2writer):
    log('Segregating reads from a merged BAM: {!r}'.format(merged_fn))

    def yield_record_and_ctg():
        """yield (r, ctg)"""
        seen = 0
        used = 0
        with io.AlignmentFile(merged_fn, 'rb', check_sq=False) as samfile:
            for r in samfile.fetch(until_eof=True):
                seen += 1
                if r.query_name not in read2ctg:
                    # print "Missing:", r.query_name
                    continue
                used += 1
                ctg = read2ctg[r.query_name]
                yield (r, ctg)
        log(' Saw {} records. Used {}.'.format(seen, used))

    # Actually write.
    for (r, ctg) in yield_record_and_ctg():
        samfn = ctg2samfn[ctg]
        #log(' Writing to samfn:{!r}'.format(samfn))
        writer = samfn2writer[samfn]
        writer.write(r)


def open_sam_writers(header, sam_fns):
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


def get_ctg2samfn(read2ctg, basedir):
    ctg2samfn = dict()
    ctgs = set(read2ctg.values())
    for ctg in ctgs:
        fn = os.path.join(basedir, ctg, '{}.bam'.format(ctg))
        ctg2samfn[ctg] = fn
    return ctg2samfn


def bam_segregate(output_fn, merged_fn):
    read2ctg_fn = merged_fn + '.read2ctg.msgpack'  # by convention
    read2ctg = io.deserialize(read2ctg_fn)
    ctg2samfn = get_ctg2samfn(read2ctg, os.path.dirname(output_fn))
    header = get_single_bam_header(merged_fn)
    bamfns = list(set(ctg2samfn.values()))
    samfn2writer = open_sam_writers(header, bamfns)
    try:
        segregate_ctgs(merged_fn, read2ctg, ctg2samfn, samfn2writer)
    finally:
        close_sam_writers(samfn2writer.values())
    fofn_content = '\n'.join(bamfns + [''])
    io.mkdirs(os.path.dirname(os.path.abspath(output_fn)))
    with open(output_fn, 'w') as ofs:
        ofs.write(fofn_content)


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='Segregate merged BAM files into single-ctg BAM files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--output-fn', type=str,
        help='A FOFN of BAM for segregated reads. The paths must encode each ctg somehow (by convention).',
    )
    parser.add_argument(
        '--merged-fn', type=str,
        help='A merged BAM file.',
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
    bam_segregate(**vars(args))


if __name__ == '__main__': # pragma: no cover
    main()
