#!/bin/bash
source bamboo_setup.sh
set -vex

python setup.py -v bdist_wheel

WHEELHOUSE=/home/cdunn/wheelhouse/gcc-6/
pip install -v --user --find-links=${WHEELHOUSE} pytest networkx pysam msgpack-python pylint
# We could try to avoid slow internet checks/downloads by using our local WHEELHOUSE, with
#   --no-index

pushd ../pypeFLOW
pip install -v --user --no-deps --edit .
popd

pushd ../FALCON
pip install -v --user --no-deps --edit .
popd

python -c 'import falcon_kit; print falcon_kit.falcon'

pip -v install --user --no-deps --use-wheel --find-links=dist/ .

pip install --user pytest pytest-cov pylint
export MY_TEST_FLAGS="-v -s --durations=0 --cov=falcon_unzip --cov-report=term-missing --cov-report=xml:coverage.xml --cov-branch"
make test
sed -i -e 's@filename="@filename="./falcon_unzip/@g' coverage.xml

pylint --errors-only falcon_unzip

ls -larth
