import argparse
import sys
from .. import dedup_h_tigs


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='Dedup h_tigs for a given ctg_id.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'ctg_id', type=str,
    )
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    args = parse_args(argv)
    dedup_h_tigs.run(**vars(args))


if __name__ == "__main__": # pragma: no cover
    main()
