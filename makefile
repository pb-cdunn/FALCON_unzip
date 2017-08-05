default:
pylint:
	pylint --errors-only falcon_unzip
test:
	py.test -v test/test_partition_and_merge_bam.py
.PHONY: test
