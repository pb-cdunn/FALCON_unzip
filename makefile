MY_TEST_FLAGS?=-vs --durations=0
WHEELHOUSE?="/mnt/software/p/python/wheelhouse/develop/"

default:
pylint:
	pylint --errors-only falcon_unzip
wheel:
	which pip
	pip wheel -v --wheel-dir=./wheelhouse --no-deps .
	ls -larth ${WHEELHOUSE}
install-wheel: wheel
	pip -v install --user --no-deps --use-wheel --find-links=dist/ .
test:
	python -c 'import falcon_kit; print falcon_kit.falcon; import falcon_unzip'
	py.test ${MY_TEST_FLAGS} --junit-xml=test.xml --doctest-modules falcon_unzip/ test/
autopep8:
	autopep8 --max-line-length=120 -ir -j0 falcon_unzip/tasks examples/ src/ setup.py test/

.PHONY: test
