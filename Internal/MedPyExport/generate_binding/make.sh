#!/bin/bash

DIST_NAME=${1-unknown}
SKIP_PROMPT=${2-0}
FOR_AM=${3-0}
if [ $DIST_NAME == "unknown" ]; then
	if [[ ${PYTHON_INCLUDE_DIR} == *"/python36/"* ]]; then DIST_NAME="medial-python36"
	elif [[ ${PYTHON_INCLUDE_DIR} == *"/python38"* ]]; then DIST_NAME="medial-python38"
	elif [[ ${PYTHON_INCLUDE_DIR} == *"/python310"* ]]; then DIST_NAME="medial-python310"
	elif [[ ${PYTHON_INCLUDE_DIR} == "/opt/medial/python27"* ]]; then DIST_NAME="medial-python27"
	elif [[ ${PYTHON_INCLUDE_DIR} == *"anaconda2"* ]]; then DIST_NAME="anaconda2"
	elif [[ ${PYTHON_INCLUDE_DIR} == "/usr"* ]]; then DIST_NAME="rh-python27"
	fi
fi

echo "(II) Python Include dir: '${PYTHON_INCLUDE_DIR}'"
echo "(II) Python Library: '${PYTHON_LIBRARY}'"
echo "(II) Compiling Python distribution: '${DIST_NAME}'"
echo "(II) AlgoMarker for client model: '${FOR_AM}'"

if [ $SKIP_PROMPT -lt 1 ]; then
	read -p "Press [Enter] to approve"
fi


if [ $FOR_AM -gt 0 ]; then
	cp SWIG.CMakeLists.AM.txt CMakeLists.txt
	cp MedPython/SWIG.CMakeLists.AM.txt MedPython/CMakeLists.txt
else
	cp SWIG.CMakeLists.txt CMakeLists.txt
	cp MedPython/SWIG.CMakeLists.txt MedPython/CMakeLists.txt
fi
mkdir -p $MR_ROOT/Libs/Internal/MedPyExport/generate_binding/CMakeBuild/Linux/Release
pushd $MR_ROOT/Libs/Internal/MedPyExport/generate_binding/CMakeBuild/Linux/Release 
cmake ../../../

version_txt=`get_git_status_text.py`
echo -e "Git version info:\n${version_txt}"
touch ${MR_ROOT}/Libs/Internal/MedUtils/MedUtils/MedGitVersion.h
if [ $FOR_AM -gt 0 ]; then
	touch ${MR_ROOT}/Libs/Internal/MedPyExport/MedPyExport/MPPidRepository.h
	make -j 20 -e GIT_HEAD_VERSION="AM_API_$version_txt";
else
	make -j 20 -e GIT_HEAD_VERSION="$version_txt";
fi
popd

NEW_RELEASE_PATH=${MR_ROOT}/Libs/Internal/MedPyExport/generate_binding/Release/${DIST_NAME}
RELEASE_PATH=${MR_ROOT}/Libs/Internal/MedPyExport/generate_binding/CMakeBuild/Linux/Release/MedPython
mkdir -p ${NEW_RELEASE_PATH}
cp ${RELEASE_PATH}/medpython.py ${RELEASE_PATH}/_medpython.so ${NEW_RELEASE_PATH}
echo "from medpython import * ; import medpython as _med ; __doc__=_med.__doc__ ; __all__=_med.__all__ ;" > ${NEW_RELEASE_PATH}/med.py
echo "Extension files copied to ${NEW_RELEASE_PATH}"
