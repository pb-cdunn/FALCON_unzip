"""
This script hacks together a JSON for the RID_TO_PHASE gathering stage.
"""
import json
import sys

def dump(sin, sout):
    data = []  # list of dict, actually

    for line in sin:
        line = line.rstrip()
        data.append(dict(rid_to_phase_out=line))
    sout.write(json.dumps(data))


def main():
    dump(sys.stdin, sys.stdout)

if __name__ == "__main__":
    main()
