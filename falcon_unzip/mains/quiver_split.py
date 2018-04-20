from falcon_kit.FastaReader import FastaReader
import argparse
import logging
import os
import sys
from .. import io
from ..tasks.quiver import TASK_QUIVER_RUN_SCRIPT

LOG = logging.getLogger()


def run(p_ctg_fasta_fn, h_ctg_fasta_fn, ctg2bamfn_fn, split_fn, bash_template_fn):
    LOG.info('Splitting quiver tasks into {!r}.'.format(
        split_fn))
    with open(bash_template_fn, 'w') as stream:
        stream.write(TASK_QUIVER_RUN_SCRIPT)
    ctg2bamfn = io.deserialize(ctg2bamfn_fn)

    ref_seq_data = {}

    # I think this will crash if the file is empty. Maybe that is ok.
    p_ctg_fa = FastaReader(p_ctg_fasta_fn)
    ctg_types = {}
    for r in p_ctg_fa:
        rid = r.name.split()[0]
        ref_seq_data[rid] = r.sequence
        ctg_types[rid] = 'p'

    # I think this will crash if the file is empty. Maybe that is ok.
    h_ctg_fa = FastaReader(h_ctg_fasta_fn)
    for r in h_ctg_fa:
        rid = r.name.split()[0]
        ref_seq_data[rid] = r.sequence
        ctg_types[rid] = 'h'

    ctg_ids = sorted(ref_seq_data.keys())
    jobs = list()
    for ctg_id in ctg_ids:
        if ctg_id not in ctg2bamfn:
            continue
        ctg_type = ctg_types[ctg_id]
        if ctg_type not in ('p', 'h'):
            msg = 'Type is {!r}, not "p" or "h". Why are we running Quiver?'.format(ctg_type)
            raise Exception(msg)
        read_bam = ctg2bamfn[ctg_id]
        if not os.path.exists(read_bam):
            # Can this ever happen?
            continue
        # The segregated *.sam were created in task_segregate.
        # Network latency should not matter (much) because we have already waited for the 'done' file.
        #m_ctg_id = ctg_id.split('-')[0] # Why did we care about this? It used to be the workdir.
        ref_fasta = os.path.abspath(os.path.join('refs', ctg_id, 'ref.fa'))
        if not os.path.exists(ref_fasta):
            io.mkdirs(os.path.dirname(ref_fasta))
            # TODO(CD): Up to 50MB of seq data. Should do this on remote host.
            #   See https://github.com/PacificBiosciences/FALCON_unzip/issues/59
            sequence = ref_seq_data[ctg_id]
            with open(ref_fasta, 'w') as f:
                print >>f, '>' + ctg_id
                print >>f, sequence
        ctg_type_fn = os.path.abspath(os.path.join('refs', ctg_id, 'ctg_type'))
        with open(ctg_type_fn, 'w') as ofs:
            ofs.write(ctg_type) # just a letter
        #wd = os.path.join(os.getcwd(), '..', 'quiver_run', ctg_id)
        cns_fasta = 'cns.fasta.gz'
        cns_fastq = 'cns.fastq.gz'
        cns_vcf = 'cns.vcf'
        job_done = 'quiver_done'
        out_ctg_type_fn = 'ctg_type'
        new_job = {}
        new_job['params'] = dict(
                ctg_id=ctg_id,
                #ctg_type=ctg_type,
        )
        new_job['input'] = dict(
                ref_fasta=ref_fasta,
                read_bam=read_bam,
                ctg_type=ctg_type_fn,
        )
        new_job['output'] = dict(
                cns_fasta=cns_fasta,
                cns_fastq=cns_fastq,
                cns_vcf=cns_vcf,
                job_done=job_done,
                ctg_type_again=out_ctg_type_fn,
        )
        new_job['wildcards'] = {
                'ctg_id': ctg_id,
                #'ctg_type': ctg_type,
        }
        #jobs[ctg_id] = new_job # TODO(CD): What about ctg_type???
        jobs.append(new_job)
    io.serialize(split_fn, jobs)


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
        '--p-ctg-fasta-fn',
        help='Typically from 3-unzip/all_p_ctg.fa',
    )
    parser.add_argument(
        '--h-ctg-fasta-fn',
        help='Typically from 3-unzip/all_h_ctg.fa',
    )
    parser.add_argument(
        '--ctg2bamfn-fn',
        help='Map of ctg -> segregated BAM (???)',
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
