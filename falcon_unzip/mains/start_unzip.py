import argparse
import sys
from .. import unzip


class HelpF(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='Run stage 3-unzip, given the results of stage 2-asm-falcon.',
        formatter_class=HelpF,
    )
    parser.add_argument(
        '--logging-config-fn',
        help='Optional standard Python logging config',
    )
    parser.add_argument(
        'config_fn', type=str,
        help='Configuration file. (This needs its own help section. Note: smrt_bin is deprecated, but if supplied will be appended to PATH.)',
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    unzip.run(**vars(args))


if __name__ == '__main__':  # pragma: no cover
    main()
