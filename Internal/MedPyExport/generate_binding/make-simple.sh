#!/bin/bash

DIST_NAME=${1-unknown}
if [ $DIST_NAME == "unknown" ]; then
	if [[ ${PYTHON_INCLUDE_DIR} == *"/python36/"* ]]; then DIST_NAME="medial-python36"
	elif [[ ${PYTHON_INCLUDE_DIR} == *"/python38"* ]]; then DIST_NAME="medial-python38"
	elif [[ ${PYTHON_INCLUDE_DIR} == *"/python310"* ]]; then 
		DIST_NAME="medial-python310"
		export CPATH=$CPATH:/nas1/Work/python-env/python310/lib/python3.10/site-packages/numpy/core/include/
		export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/nas1/Work/python-env/python310/lib/python3.10/site-packages/numpy/core/include/
	else
		echo "Unknowwn python version"
		exit -1
	fi
fi

echo "(II) Python Include dir: '${PYTHON_INCLUDE_DIR}'"
echo "(II) Python Library: '${PYTHON_LIBRARY}'"
echo "(II) Compiling Python distribution: '${DIST_NAME}'"

COMPILE_STATIC_PATH=${MR_ROOT}/Libs/Internal/AlgoMarker/LinuxSharedLib/Build.ubuntu/lib
#STATIC_LIBS_TARGET=/nas1/Work/SharedLibs/linux/ubuntu/static_libs/Release
STATIC_LIBS_TARGET=/nas1/Work/SharedLibs/linux/ubuntu/static_libs/Release_new

if [[ "`which sudo 2>>/dev/null`" ]]; 
then 
export SUDO="`which sudo`"
else
export SUDO=''
fi

mkdir -p `realpath ${COMPILE_STATIC_PATH}/..`
mkdir -p ${COMPILE_STATIC_PATH}
echo "STATIC_LIBS_TARGET=${STATIC_LIBS_TARGET}"
echo "COMPILE_STATIC_PATH=${COMPILE_STATIC_PATH}"

cp -R ${STATIC_LIBS_TARGET}/* ${COMPILE_STATIC_PATH}

source /nas1/Work/python-env/python310/bin/activate

cp SWIG.CMakeLists.txt CMakeLists.txt
cp MedPython/SWIG.CMakeLists.txt MedPython/CMakeLists.txt

mkdir -p $MR_ROOT/Libs/Internal/MedPyExport/generate_binding/CMakeBuild/Linux/Release
pushd $MR_ROOT/Libs/Internal/MedPyExport/generate_binding/CMakeBuild/Linux/Release 
cmake ../../../

version_txt=`get_git_status_text.py`
echo -e "Git version info:\n${version_txt}"
touch ${MR_ROOT}/Libs/Internal/MedUtils/MedUtils/MedGitVersion.h

make -j 20 -e GIT_HEAD_VERSION="$version_txt";

popd

NEW_RELEASE_PATH=${MR_ROOT}/Libs/Internal/MedPyExport/generate_binding/Release/${DIST_NAME}
RELEASE_PATH=${MR_ROOT}/Libs/Internal/MedPyExport/generate_binding/CMakeBuild/Linux/Release/MedPython
mkdir -p ${NEW_RELEASE_PATH}
cp ${RELEASE_PATH}/medpython.py ${RELEASE_PATH}/_medpython.so ${NEW_RELEASE_PATH}
echo "from medpython import * ; import medpython as _med ; __doc__=_med.__doc__ ; __all__=_med.__all__ ;" > ${NEW_RELEASE_PATH}/med.py
echo "Extension files copied to ${NEW_RELEASE_PATH}"
