"""NOT CURRENTLY USED
"""
from ..io import (serialize, yield_bam_fn, log, AlignmentFile)
import os


def select_reads_from_bam(input_bam_fofn_fn, rawread_to_contigs_fn, rawread_ids_fn, sam_dir,
                          max_n_open_files):
    """Write ctg.sam files into cwd,
    for each 'ctg' read in input BAMs.
    """
    read_partition = {}
    read_to_ctgs = {}

    log("rawread_ids_fn:", repr(rawread_ids_fn))
    log("rawread_to_contigs_fn:", repr(rawread_to_contigs_fn))
    rid_to_oid = open(rawread_ids_fn).read().split('\n')
    with open(rawread_to_contigs_fn) as f:
        for row in f:
            row = row.strip().split()
            if int(row[3]) >= 1:  # keep top one hits
                continue
            ctg_id = row[1]
            if ctg_id == 'NA':
                continue
            read_partition.setdefault(ctg_id, set())
            r_id = row[0]
            o_id = rid_to_oid[int(r_id)]
            read_partition[ctg_id].add(o_id)
            read_to_ctgs.setdefault(o_id, [])
            read_to_ctgs[o_id].append((int(row[4]), ctg_id))
    log("num read_partitions:", len(read_partition))
    log("num read_to_ctgs:", len(read_to_ctgs))

    header = None
    for fn in yield_bam_fn(input_bam_fofn_fn):
        with AlignmentFile(fn, 'rb', check_sq=False) as samfile:
            if header == None:
                header = samfile.header
            else:
                header['RG'].extend(samfile.header['RG'])
    try:
        PG = header.pop('PG')  # remove PG line as there might be a bug that generates no readable chrs
    except KeyError:
        pass
    # print PG

    #base_dir = os.getcwd()
    #outfile = AlignmentFile( os.path.join(base_dir, 'header.sam' ), 'wh', header=header )
    # outfile.close()

    ctgs = read_partition.keys()
    ctgs.sort()
    selected_ctgs = set()
    for ctg in ctgs:
        picked_reads = read_partition[ctg]
        # print "ctg, len:", ctg, len(picked_reads)
        if len(picked_reads) > 20:
            selected_ctgs.add(ctg)

    merge_and_split_alignments(input_bam_fofn_fn, read_to_ctgs, selected_ctgs, header, sam_dir,
                               max_n_open_samfiles=max(1, max_n_open_files - 1))


def merge_and_split_alignments(input_bam_fofn_fn, read_to_ctgs, selected_ctgs, header, sam_dir, max_n_open_samfiles):
    """
    For each AlignmentFile in input_bam_fofn_fn,
      for each record in that file,
        find the ctgs from the query_name (based on read_to_ctgs and selected_ctgs),
        choose the first ctg (lexigraphically),
        and append the record to a new samfile (in sam_dir) based on the ctg name.

    That is the logic. In reality, to avoid too many open samfiles, we may use extra passes.
    When max_n_open_samfiles==1, this is O(n*m), for n ctgs and m BAM alignment records; we scan
    each BAM input n/max_n_open_samfiles times.
    """
    def yield_record_and_ctg():
        """yield (r, ctg)"""
        for fn in yield_bam_fn(input_bam_fofn_fn):
            with AlignmentFile(fn, 'rb', check_sq=False) as samfile:
                for r in samfile.fetch(until_eof=True):
                    if r.query_name not in read_to_ctgs:
                        # print "Missing:", r.query_name
                        continue
                    ctg_list = read_to_ctgs[r.query_name]
                    ctg_list.sort()
                    score, ctg = ctg_list[0]
                    if ctg not in selected_ctgs:
                        # print "Not selected:", ctg
                        continue
                    yield r, ctg

    outfilenames = {}
    ctgs = list()

    for (r, ctg) in yield_record_and_ctg():
        if ctg not in outfilenames:
            samfile_fn = os.path.join(sam_dir, '%s.bam' % ctg)
            log('samfile_fn:{!r}'.format(samfile_fn))
            outfilenames[ctg] = samfile_fn
            ctgs.append(ctg)
    ctgs.sort(reverse=True)  # Order does not matter, but for debugging we might as well sort.

    while ctgs:
        ctg_openset = set()
        while ctgs and len(ctg_openset) < max_n_open_samfiles:
            ctg_openset.add(ctgs.pop())
        log('ctg_openset:', ctg_openset)
        outfile = {}
        for (r, ctg) in yield_record_and_ctg():
            if ctg not in ctg_openset:
                continue
            samfile_fn = outfilenames[ctg]
            if ctg not in outfile:
                log('Opening samfile_fn:{!r}'.format(samfile_fn))
                outfile[ctg] = AlignmentFile(samfile_fn, 'wb', header=header)
            #log('Writing to samfile_fn:{!r}'.format(samfile_fn))
            outfile[ctg].write(r)
        for ctg in outfile:
            outfile[ctg].close()


##################################
import argparse
import logging
import sys


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='Write ctg.sam files, based on BAM subreads.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--rawread-to-contigs', type=str,
        default='./2-asm-falcon/read_maps/dump_rawread_ids/rawread_to_contigs', help='rawread_to_contigs file (from where?)')
    parser.add_argument(
        '--rawread-ids', type=str,
        default='./2-asm-falcon/read_maps/dump_rawread_ids/rawread_ids', help='rawread_ids file (from where?)')
    parser.add_argument(
        '--sam-dir', type=str, default='./4-quiver/reads', help='Output directory for ctg.sam files')
    parser.add_argument(
        '--max-n-open-files', type=int, default='50',
        help='We write sam files several at-a-time, limited by this.')
    parser.add_argument(
        'input_bam_fofn', type=str,
        help='File of BAM filenames. Paths are relative to dir of FOFN, not CWD.')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(message)s',
    )

    select_reads_from_bam(args.input_bam_fofn, args.rawread_to_contigs, args.rawread_ids, args.sam_dir,
                          args.max_n_open_files)


if __name__ == '__main__':  # pragma: no cover
    main()
