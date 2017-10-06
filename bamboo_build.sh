#!/bin/bash
source bamboo_setup.sh
set -vex

python setup.py -v bdist_wheel

WHEELHOUSE=/home/cdunn/wheelhouse/gcc-6/
pip install --user --find-links=${WHEELHOUSE} pytest networkx pysam msgpack-python pylint

pushd ../pypeFLOW
pip install -v --user --no-deps --edit .
popd

pushd ../FALCON
pip install -v --user --no-deps --edit .
popd

pip -v install --user --no-deps --use-wheel --find-links=dist/ .

mkdir -p test-reports/
export MY_TEST_FLAGS

MY_TEST_FLAGS="-vs --durations=0 --junitxml=test-reports/pytest.xml"
py.test ${MY_TEST_FLAGS} test/

MY_TEST_FLAGS="-vs --durations=0 --junitxml=test-reports/doctest.xml"
py.test ${MY_TEST_FLAGS} --doctest-modules falcon_unzip/

pylint --errors-only falcon_unzip
