default:
pylint:
	pylint --errors-only falcon_unzip
test:
	python -c 'import falcon_kit; print falcon_kit.falcon'
	#pip install --user pytest
	py.test ${MY_TEST_FLAGS} test/

.PHONY: test
