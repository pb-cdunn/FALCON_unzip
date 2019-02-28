import argparse
import logging
import os
import sys
from .. import io
from ..tasks.unzip import TASK_SEGREGATE_RUN_SCRIPT

LOG = logging.getLogger()


def run(unzip_config_fn, merged_fofn_fn, bash_template_fn, split_fn):
    LOG.info('Scatting segregate-reads tasks into {!r}. Reading merged BAM names from FOFN: {!r}'.format(
        split_fn, merged_fofn_fn))
    with open(bash_template_fn, 'w') as stream:
        stream.write(TASK_SEGREGATE_RUN_SCRIPT)

    unzip_config = io.deserialize(unzip_config_fn)
    extra = '1' if unzip_config['polish_include_zmw_all_subreads'] else '0'

    basedir = os.path.normpath(os.path.dirname(split_fn))
    fns = list(io.yield_abspath_from_fofn(merged_fofn_fn))
    # ctg is encoded into each filepath within each FOFN, fwiw.
    jobs = list()
    for i, merged_bam_fn in enumerate(fns):
        job_name = 'segr{:03d}'.format(i) # wildcard value # TODO: Needed?

        #segregated_bam_fofn = '{}/{}/segregated_bam.fofn'.format(
        #        basedir, job_name)
        segregated_bam_fns_fn = 'segregated-bam-fns.json'
        job = dict()
        job['input'] = dict(
                # The other input is next to this one, postfixed by convention: 'merged.bam.read2ctg.json'
                merged_bam_fn=merged_bam_fn,
        )
        job['output'] = dict(
                segregated_bam_fns=segregated_bam_fns_fn,
        )
        job['params'] = dict(
                extra=extra,
        )
        job['wildcards'] = {'segr': job_name} # This should match the wildcard used in the pattern elsewhere.
        jobs.append(job)
    io.serialize(split_fn, jobs)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Split the tasks which generate segregated BAM.'
    epilog = 'To learn about inputs and outputs in the serialized jobs, see bam_segregate.py'
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--unzip-config-fn', required=True,
        help='Input. Serialized [Unzip] section of config.',
    )
    parser.add_argument(
        '--merged-fofn-fn',
        help='Input. FOFN of merged BAM.',
    )
    parser.add_argument(
        '--split-fn',
        help='Output. JSON list of units of work.')
    parser.add_argument(
        '--bash-template-fn',
        help='Output. Copy of known bash template, for use later.')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
