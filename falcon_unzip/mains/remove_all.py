from .. import io
import argparse
import logging
import os
import sys

LOG = logging.getLogger(__name__)


def run(fofn_fn, fail_on_error):
    fofn = io.deserialize(fofn_fn)
    for fn in fofn:
        try:
            os.unlink(fn)
        except Exception:
            if fail_on_error:
                raise
            LOG.warning('Could not rm {!r}. Ignoring.'.format(fn))


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='Remove all files in a serialized list of filenames.',
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--fofn-fn', required=True,
        help='Input. JSON, msgpack, or FOFN.')
    parser.add_argument(
        '--fail-on-error', action='store_true',
        help='Normally, we skip errors, since at worst we leave a file sitting there.')
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    logging.basicConfig(level=logging.INFO)
    main()
