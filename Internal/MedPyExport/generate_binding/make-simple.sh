#!/bin/bash
set -e
DIST_NAME=${1-unknown}
PY_VERSION=$(python --version | awk '{print $2}' | awk -F. '{print $1 "." $2}')
PY_VERSION_SHORT=$(python --version | awk '{print $2}' | awk -F. '{print $1 $2}')
CURRENT_DIR=$(realpath ${0%/*})

if [ $DIST_NAME == "unknown" ]; then
	DIST_NAME="medial-python${PY_VERSION_SHORT}"
fi

#echo "(II) Python Include dir: '${PYTHON_INCLUDE_DIR}'"
#echo "(II) Python Library: '${PYTHON_LIBRARY}'"
echo "(II) Compiling Python distribution: '${DIST_NAME}'"

#STATIC_LIBS_TARGET=/nas1/Work/SharedLibs/linux/ubuntu/static_libs/Release_new

mkdir -p ${CURRENT_DIR}/CMakeBuild/Linux/Release
pushd ${CURRENT_DIR}/CMakeBuild/Linux/Release 
cmake ../../../

set +e
GIT_COMMIT_HASH=$(git rev-parse HEAD)
version_txt=$(date +'Commit_'${GIT_COMMIT_HASH}'_Build_On_%Y%m%d_%H:%M:%S')
set -e
echo -e "Git version info:\n${version_txt}"
touch ${CURRENT_DIR}/../../MedUtils/MedUtils/MedGitVersion.h

make -j $(nproc) -e GIT_HEAD_VERSION="$version_txt";

popd

NEW_RELEASE_PATH=${CURRENT_DIR}/Release/${DIST_NAME}
RELEASE_PATH=${CURRENT_DIR}/CMakeBuild/Linux/Release/MedPython
mkdir -p ${NEW_RELEASE_PATH}
cp ${RELEASE_PATH}/../medpython.py ${RELEASE_PATH}/_medpython.so ${NEW_RELEASE_PATH}
echo "from medpython import * ; import medpython as _med ; __doc__=_med.__doc__ ; __all__=_med.__all__ ;" > ${NEW_RELEASE_PATH}/med.py
echo "Extension files copied to ${NEW_RELEASE_PATH}"