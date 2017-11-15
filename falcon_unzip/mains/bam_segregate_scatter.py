import argparse
import logging
import os
import sys
from .. import io

LOG = logging.getLogger()


def run(merged_fofn_fn, scattered_fn):
    LOG.info('Scatting segregate-reads tasks into {!r}. Reading merged BAM names from FOFN: {!r}'.format(
        scattered_fn, merged_fofn_fn))
    basedir = os.path.dirname(os.path.abspath(scattered_fn))
    fns = list(io.yield_abspath_from_fofn(merged_fofn_fn))
    # ctg is encoded into each filepath within each FOFN, fwiw.
    jobs = list()
    for i, merged_bamfn in enumerate(fns):
        job_name = 'segr{:03d}'.format(i) # wildcard value

        segregated_bam_fofn = '{}/{}/segregated_bam.fofn'.format(
                basedir, job_name)
        job = dict()
        job['input'] = dict(
                # The other input is next to this one, postfixed by convention: 'merged.bam.read2ctg.json'
                merged_bamfn=merged_bamfn,
        )
        job['output'] = dict(
                segregated_bam_fofn=segregated_bam_fofn,
        )
        job['params'] = dict()
        job['wildcards'] = {'segr': job_name} # This should match the wildcard used in the pattern elsewhere.
        jobs.append(job)
    io.serialize(scattered_fn, jobs)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Scatter the tasks which generate segregated BAM.'
    epilog = 'To learn about inputs and outputs in the serialized jobs, see bam_segregate.py'
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--merged-fofn-fn',
        help='FOFN of merged ???',
    )
    parser.add_argument(
        '--scattered-fn',
        help='JSON file of the ??? files for segregation',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
