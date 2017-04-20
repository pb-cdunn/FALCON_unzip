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
#job_type = SGE
#job_type = local
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
smrt_bin=/pbi/dept/secondary/builds/4.1.0/current_smrttools_incremental_installdir/smrtcmds/bin
sge_phasing= -pe smp 12 -q bigmem
sge_quiver= -pe smp 12 -q sequel-farm
sge_track_reads= -pe smp 12 -q default
sge_blasr_aln= -pe smp 24 -q bigmem
sge_hasm= -pe smp 48 -q bigmem
unzip_concurrent_jobs = 12
quiver_concurrent_jobs = 12
```
