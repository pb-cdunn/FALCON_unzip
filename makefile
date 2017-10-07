MY_TEST_FLAGS?=-vs --durations=0

default:
pylint:
	pylint --errors-only falcon_unzip
test:
	python -c 'import falcon_kit; print falcon_kit.falcon'
	py.test ${MY_TEST_FLAGS} test/
doctest:
	py.test ${MY_TEST_FLAGS} --doctest-modules falcon_unzip/
wheel:
	# This is basically a syntax check.
	python setup.py bdist_wheel
install-wheel: wheel
	pip -v install --user --no-deps --use-wheel --find-links=dist/ .

.PHONY: test
