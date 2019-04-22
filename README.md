Please refer to https://github.com/PacificBiosciences/FALCON_unzip/wiki

Needs (in `smrt_bin` directory):

* blasr
* samtools
* nucmer
* variantCaller.py

Tested with smrttools-4.1.0

## Running
For an example:
* https://github.com/pb-cdunn/FALCON-examples/tree/master/run/greg200k-sv2

But basically, this works:
```sh
# after running fc_run.py from FALCON,
fc_unzip.py fc_unzip.cfg
fc_quiver.py fc_unzip.cfg
```

## Example config

These are the available settings:
```ini
[General]
#use_tmpdir = /scratch
pwatcher_type = blocking
job_type = string
job_queue = bash -C ${CMD} >| ${STDOUT_FILE} 2>| ${STDERR_FILE}

#job_queue = bash -C ${CMD}
# By dropping STD*_FILE, we see all output on the console.
# That helps debugging in TravisCI/Bamboo.

max_n_open_files = 1000

[Unzip]
input_fofn= input.fofn
input_bam_fofn= input_bam.fofn
#polish_vc_ignore_error = true
#polish_use_blasr = true
#polish_include_zmw_all_subreads = true

[job.step.unzip.track_reads]
njobs=1
NPROC=48

[job.step.unzip.blasr_aln]
njobs=8
NPROC=16

[job.step.unzip.phasing]
njobs=16
NPROC=2

[job.step.unzip.hasm]
njobs=1
NPROC=48

[job.step.unzip.quiver]
njobs=12
NPROC=12
```
For job-distribution, see https://github.com/PacificBiosciences/pypeFLOW/wiki/Configuration
