import argparse
import logging
import os
import sys
from .. import io

LOG = logging.getLogger()


def run(gathered_fn, rid_to_phase_all_fn):
    gathered_dir = os.path.dirname(gathered_fn)
    out_dir = os.path.dirname(rid_to_phase_all_fn)
    rid_to_phase_all = list()
    gathered = io.deserialize(gathered_fn)
    for job_output in gathered.values():
        fn = job_output['fns']['rid_to_phase_out']
        if not os.path.isabs(fn):
            fn = os.path.relpath(os.path.join(gathered_dir, fn), out_dir)
        # Discern ctgs from filepaths. (Or we could we use the key in gathered.json.)
        ctg_id = os.path.basename(os.path.dirname(fn))
        with open(fn) as stream:
            for line in stream:
                rid_to_phase_all.append(line.strip())
    with open(rid_to_phase_all_fn, 'w') as stream:
        stream.write('\n'.join(rid_to_phase_all))
        stream.write('\n')


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    description = 'Concatenate rid_to_phase files.'
    epilog = ''
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--gathered-fn',
        help='Input: dict of ctg-> dict of "fn"->rid_to_phase filename',
    )
    parser.add_argument(
        '--rid-to-phase-all-fn',
        help='Output: Concatenation of rid_to_phase files',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
