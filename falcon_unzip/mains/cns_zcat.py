"""TODO: This can be a script, after we figure out how to deal with scatter/gather.
"""
import argparse
import logging
import sys
from .. import io


def run(rm_intermediates,
        gathered_quiver_fn,
        cns_p_ctg_fasta_fn,
        cns_p_ctg_fastq_fn,
        cns_h_ctg_fasta_fn,
        cns_h_ctg_fastq_fn,
):
    gathered_dict = io.deserialize(gathered_quiver_fn)
    # 'p_ctg': list(sorted(p_ctg_out)), # cns_fasta_fn, cns_fastq_fn
    # 'h_ctg': list(sorted(h_ctg_out)), # cns_fasta_fn, cns_fastq_fn

    used_files = list()

    io.rm(cns_p_ctg_fasta_fn)
    io.touch(cns_p_ctg_fasta_fn)
    io.rm(cns_p_ctg_fastq_fn)
    io.touch(cns_p_ctg_fastq_fn)
    for cns_fasta_fn, cns_fastq_fn in gathered_dict['p_ctg']:
            io.syscall('zcat {cns_fasta_fn} >> {cns_p_ctg_fasta_fn}'.format(**locals()))
            io.syscall('zcat {cns_fastq_fn} >> {cns_p_ctg_fastq_fn}'.format(**locals()))
            used_files.append(cns_fasta_fn)
            used_files.append(cns_fastq_fn)

    io.rm(cns_h_ctg_fasta_fn)
    io.touch(cns_h_ctg_fasta_fn)
    io.rm(cns_h_ctg_fastq_fn)
    io.touch(cns_h_ctg_fastq_fn)
    for cns_fasta_fn, cns_fastq_fn in gathered_dict['h_ctg']:
            io.syscall('zcat {cns_fasta_fn} >> {cns_h_ctg_fasta_fn}'.format(**locals()))
            io.syscall('zcat {cns_fastq_fn} >> {cns_h_ctg_fastq_fn}'.format(**locals()))

    if rm_intermediates:
        io.rm(*used_files)

class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Concatenate results of consensus, and zip, then (maybe) rm intermediates'
    epilog = 'The FOFNs must have a fasta and a fastq filename on each line.'
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--rm-intermediates', action='store_true',
    )
    parser.add_argument(
        '--gathered-quiver-fn',
        required=True,
        help='A very specific format. See yield_quiver_tasks() in this repo.',
    )
    parser.add_argument(
        '--cns-p-ctg-fasta-fn',
        required=True,
        help='output',
    )
    parser.add_argument(
        '--cns-p-ctg-fastq-fn',
        required=True,
        help='output',
    )
    parser.add_argument(
        '--cns-h-ctg-fasta-fn',
        required=True,
        help='output',
    )
    parser.add_argument(
        '--cns-h-ctg-fastq-fn',
        required=True,
        help='output',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
