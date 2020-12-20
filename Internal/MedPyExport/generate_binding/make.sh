#!/bin/bash

DIST_NAME=${1-unknown}
if [ $DIST_NAME == "unknown" ]; then
	if [[ ${PYTHON_INCLUDE_DIR} == *"/python36/"* ]]; then DIST_NAME="medial-python36"
	elif [[ ${PYTHON_INCLUDE_DIR} == "/opt/medial/python27"* ]]; then DIST_NAME="medial-python27"
	elif [[ ${PYTHON_INCLUDE_DIR} == *"anaconda2"* ]]; then DIST_NAME="anaconda2"
	elif [[ ${PYTHON_INCLUDE_DIR} == "/usr"* ]]; then DIST_NAME="rh-python27"
	fi
fi

echo "(II) Python Include dir: '${PYTHON_INCLUDE_DIR}'"
echo "(II) Python Library: '${PYTHON_LIBRARY}'"
echo "(II) Compiling Python distribution: '${DIST_NAME}'"

read -p "Press [Enter] to approve"

cp SWIG.CMakeLists.txt CMakeLists.txt
cp MedPython/SWIG.CMakeLists.txt MedPython/CMakeLists.txt
mkdir -p $MR_ROOT/Libs/Internal/MedPyExport/generate_binding/CMakeBuild/Linux/Release
pushd $MR_ROOT/Libs/Internal/MedPyExport/generate_binding/CMakeBuild/Linux/Release 
cmake ../../../
make -j 8;
popd

NEW_RELEASE_PATH=${MR_ROOT}/Libs/Internal/MedPyExport/generate_binding/Release/${DIST_NAME}
RELEASE_PATH=${MR_ROOT}/Libs/Internal/MedPyExport/generate_binding/CMakeBuild/Linux/Release/MedPython
mkdir -p ${NEW_RELEASE_PATH}
cp ${RELEASE_PATH}/medpython.py ${RELEASE_PATH}/_medpython.so ${NEW_RELEASE_PATH}
echo "from medpython import * ; import medpython as _med ; __doc__=_med.__doc__ ; __all__=_med.__all__ ;" > ${NEW_RELEASE_PATH}/med.py
echo "Extension files copied to ${NEW_RELEASE_PATH}"
