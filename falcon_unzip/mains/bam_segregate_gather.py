import argparse
import logging
import os
import sys
from .. import io

LOG = logging.getLogger()


def run(gathered_fn, ctg2segregated_bamfn_fn):
    ctg2segregated_bamfn = dict()
    gathered = io.deserialize(gathered_fn)
    for job_output in gathered.values():
        # Read FOFN.
        fofn_fn = job_output['fns']['segregated_bam_fofn']
        segregated_bam_fns = list(io.yield_abspath_from_fofn(fofn_fn))
        # Discern ctgs from filepaths.
        for bamfn in segregated_bam_fns:
            basename = os.path.basename(bamfn)
            ctg = os.path.splitext(basename)[0]
            ctg2segregated_bamfn[ctg] = bamfn
    io.serialize(ctg2segregated_bamfn_fn, ctg2segregated_bamfn)
    io.serialize(ctg2segregated_bamfn_fn + '.json', ctg2segregated_bamfn)  # for debugging


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Gather segregated BAM files'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--gathered-fn',
        help='Input: serialized something (list of BAM-fofn filenames?)',
    )
    parser.add_argument(
        '--ctg2segregated-bamfn-fn',
        help='Output: serialized map of ctg -> BAM-fn',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
