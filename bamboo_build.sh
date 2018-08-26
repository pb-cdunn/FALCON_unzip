#!/bin/bash
source bamboo_setup.sh
set -vex

python setup.py -v bdist_wheel

WHEELHOUSE="/mnt/software/p/python/wheelhouse/develop/"
pip install -v --user --no-index --find-links=${WHEELHOUSE} pytest networkx pysam msgpack pylint future intervaltree pypeflow falcon_kit

#pushd ../pypeFLOW
#pip install -v --user --no-deps --edit .
#popd
#
#pushd ../FALCON
#pip install -v --user --no-deps --edit .
#popd

python -c 'import falcon_kit; print falcon_kit.falcon'

pip -v install --user --no-deps --use-wheel --find-links=${WHEELHOUSE} .

pip install --user pytest pytest-cov pylint
export MY_TEST_FLAGS="-v -s --durations=0 --cov=falcon_unzip --cov-report=term-missing --cov-report=xml:coverage.xml --cov-branch"
make test
#sed -i -e 's@filename="@filename="./falcon_unzip/@g' coverage.xml

pylint --errors-only falcon_unzip

ls -larth

bash bamboo_wheel.sh
