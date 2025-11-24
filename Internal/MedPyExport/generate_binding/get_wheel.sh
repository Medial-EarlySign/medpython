#!/bin/bash
CURRENT_DIR=$(realpath ${0%/*})

python -m pip install build auditwheel
mkdir -p wheelhouse
#rm -fr build dist

set +e
GIT_COMMIT_HASH=$(git rev-parse HEAD)
version_txt=$(date +'Commit_'${GIT_COMMIT_HASH}'_Build_On_%Y%m%d_%H:%M:%S')
set -e
echo -e "Git version info:\n${version_txt}"
touch ${CURRENT_DIR}/../../MedUtils/MedUtils/MedGitVersion.h
export GIT_HEAD_VERSION=$version_txt 

python -m build --wheel --outdir dist
#python -m build --sdist --outdir wheelhouse/

for whl in dist/*.whl; do
    auditwheel repair "$whl" --plat manylinux2014_x86_64 -w wheelhouse/
done