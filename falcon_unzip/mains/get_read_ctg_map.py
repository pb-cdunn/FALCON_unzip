"""Obsolete.
The tasks are now part of unzip.
"""
from __future__ import absolute_import
from ..tasks import unzip as tasks_unzip
# pylint: disable=no-name-in-module, import-error, fixme, line-too-long
from pypeflow.simple_pwatcher_bridge import (PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase,
                                             makePypeLocalFile, fn, PypeTask)
PypeThreadTaskBase = MyFakePypeThreadTaskBase
import argparse
import glob
import logging
import sys
import subprocess as sp
import shlex
import os

LOG = logging.getLogger(__name__)


def get_read_ctg_map():
    wf = PypeProcWatcherWorkflow(
        max_jobs=12,
    )
    """
            job_type=config['job_type'],
            job_queue=config['job_queue'],
            sge_option=config.get('sge_option', ''),
            watcher_type=config['pwatcher_type'],
            watcher_directory=config['pwatcher_directory'])
    """
    read_to_contig_map_file = makePypeLocalFile('3-unzip/reads/get_read_ctg_map/read_to_contig_map')
    wf.addTasks(tasks_unzip.create_tasks_read_to_contig_map(read_to_contig_map_file))
    wf.refreshTargets()  # block


def parse_args(argv):
    description = """
Based on multiple Falcon inputs,
generate "./3-unzip/reads/get_read_ctg_map/read_to_contig_map" that contains the
information from the chain of mapping: (contig id) -> (internal p-read id) -> (internal raw-read id) -> (original read id)
"""
    epilog = 'Run in falcon base-dir. Task-dirs are by convention.'
    parser = argparse.ArgumentParser(description=description,
                                     epilog=epilog,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    args = parser.parse_args(argv[1:])
    return args


def main(argv=sys.argv):
    logging.basicConfig()
    args = parse_args(argv)
    get_read_ctg_map()


if __name__ == '__main__':
    main()
