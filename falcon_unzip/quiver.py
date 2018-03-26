from falcon_kit import run_support as support
from falcon_kit import snakemake
from pypeflow.simple_pwatcher_bridge import (
    PypeProcWatcherWorkflow)
from .tasks import quiver as tasks_quiver
from . import io
from . import unzip
import glob
import logging
import os
import re
import time
import ConfigParser

LOG = logging.getLogger(__name__)


def run(config_fn, logging_config_fn):
    global LOG
    LOG = support.setup_logger(logging_config_fn)

    config = unzip.parse_config(config_fn)

    #config_absbasedir = os.path.dirname(os.path.abspath(config_fn))
    #input_bam_fofn = 'input_bam.fofn'
    #if config.has_option('Unzip', 'input_bam_fofn'):
    #    input_bam_fofn = config.get('Unzip', 'input_bam_fofn')
    #if not os.path.isabs(input_bam_fofn):
    #    input_bam_fofn = os.path.join(config_absbasedir, input_bam_fofn)

    wf = PypeProcWatcherWorkflow(
        job_defaults=config['job.defaults'],
    )
    with open('foo.snake', 'w') as snakemake_writer:
        rule_writer = snakemake.SnakemakeRuleWriter(snakemake_writer)
        tasks_quiver.run_workflow(wf, config, rule_writer)
