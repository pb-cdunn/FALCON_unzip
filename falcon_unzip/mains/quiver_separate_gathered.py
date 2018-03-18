import argparse
import logging
import os
import sys
from .. import io

LOG = logging.getLogger()


def run(gathered_fn, output_fn):
    """We wrote the "gathered" files during task construction.
    """
    p_ctg_out = list()
    h_ctg_out = list()
    gathered = io.deserialize(gathered_fn)

    def resolved(fn):
        if os.path.isabs(fn):
            # We should not need this ever, but for now we are flexible.
            return fn
        return os.path.join(os.path.dirname(gathered_fn), fn)

    for job_output in gathered:
        cns_fasta_fn = resolved(job_output['cns_fasta'])
        cns_fastq_fn = resolved(job_output['cns_fastq'])
        ctg_type_fn = resolved(job_output['ctg_type_again'])
        ctg_type = open(ctg_type_fn).read().strip()
        assert ctg_type in 'ph', 'ctg_type={!r}'.format(ctg_type)
        # I don't know what this was supposed to do. Delete it later.
        #wd = os.path.dirname(cns_fasta_fn)
        #calc_cns_fasta = '{wd}/cns.fasta.gz'.format(**locals())
        #calc_cns_fastq = '{wd}/cns.fastq.gz'.format(**locals())
        #assert cns_fasta_fn == calc_cns_fasta, '\n{!r} !=\n{!r}'.format(
        #        cns_fasta_fn, calc_cns_fasta)
        if ctg_type == 'p':
            p_ctg_out.append([cns_fasta_fn, cns_fastq_fn])
        elif ctg_type == 'h':
            h_ctg_out.append([cns_fasta_fn, cns_fastq_fn])
    output_dict = {
            'p_ctg': list(sorted(p_ctg_out)), # cns_fasta_fn, cns_fastq_fn
            'h_ctg': list(sorted(h_ctg_out)), # cns_fasta_fn, cns_fastq_fn
    }
    io.serialize(output_fn, output_dict)


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Scatter to prepare for parallel Quiver jobs'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--gathered-fn',
        help='Result of quiver-gather',
    )
    parser.add_argument(
        '--output-fn',
        help='Serialized (JSON or msgpack) file of {"p_ctg"|"h_ctg": [fasta, fastq]}',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
