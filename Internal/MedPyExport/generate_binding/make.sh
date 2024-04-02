#!/bin/bash

DIST_NAME=${1-unknown}
SKIP_PROMPT=${2-0}
if [ $DIST_NAME == "unknown" ]; then
	if [[ ${PYTHON_INCLUDE_DIR} == *"/python36/"* ]]; then DIST_NAME="medial-python36"
	elif [[ ${PYTHON_INCLUDE_DIR} == *"/python38"* ]]; then DIST_NAME="medial-python38"
	elif [[ ${PYTHON_INCLUDE_DIR} == *"/python310"* ]]; then 
		DIST_NAME="medial-python310"
		export CPATH=$CPATH:/nas1/Work/python-env/python310/lib/python3.10/site-packages/numpy/core/include/
		export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/nas1/Work/python-env/python310/lib/python3.10/site-packages/numpy/core/include/
	elif [[ ${PYTHON_INCLUDE_DIR} == "/opt/medial/python27"* ]]; then DIST_NAME="medial-python27"
	elif [[ ${PYTHON_INCLUDE_DIR} == *"anaconda2"* ]]; then DIST_NAME="anaconda2"
	elif [[ ${PYTHON_INCLUDE_DIR} == "/usr"* ]]; then DIST_NAME="rh-python27"
	fi
fi

echo "(II) Python Include dir: '${PYTHON_INCLUDE_DIR}'"
echo "(II) Python Library: '${PYTHON_LIBRARY}'"
echo "(II) Compiling Python distribution: '${DIST_NAME}'"

#NEED To Build AlgoMarker before to have libraries

if [ $SKIP_PROMPT -lt 1 ]; then
	read -p "Press [Enter] to approve"
fi

BUILD_BOOST=true
BUILD_XGBOOST=true
BUILD_LIGHTGBM=true
BUILD_ALGOMARKER=true

OVERRIDE_ALL=${3-0}

BOOST_FPIC_DIR=/nas1/Work/Libs/Boost/boost_1_67_0-fPIC.ubuntu
DIR=$MR_ROOT/Libs/Internal/AlgoMarker/LinuxSharedLib
BLDDIR="${DIR}/Build.ubuntu"
STATIC_LIBS_TARGET="${BLDDIR}/lib"

mkdir -p ${BLDDIR}
mkdir -p ${STATIC_LIBS_TARGET}

if [[ "`which sudo 2>>/dev/null`" ]]; 
then 
export SUDO="`which sudo`"
else
export SUDO=''
fi

echo "SUDO=${SUDO}"
echo "BOOST_FPIC_DIR=${BOOST_FPIC_DIR}"
echo "DIR=${DIR}"
echo "BLDDIR=${BLDDIR}"
echo "STATIC_LIBS_TARGET=${STATIC_LIBS_TARGET}"

mkdir -p "${STATIC_LIBS_TARGET}"
mkdir -p ${BOOST_FPIC_DIR}/installation
source /nas1/Work/python-env/python310/bin/activate

if [ "$BUILD_BOOST" = true ] ; then
pushd ${BOOST_FPIC_DIR}
if [ ! -f ${STATIC_LIBS_TARGET}/libboost_system.a ] || [ $OVERRIDE_ALL -gt 0 ]; then
	./b2 --clean
	./bootstrap.sh --without-icu --with-libraries=filesystem,system,regex,program_options --prefix=${BOOST_FPIC_DIR}/installation
	./b2 link=static variant=release linkflags=-static-libstdc++ -j20 cxxflags="-fPIC" toolset=gcc --stagedir="${BLDDIR}"
else
	echo "BOOST Already exists"
fi
popd
fi

if [ "$BUILD_XGBOOST" = true ] ; then
	pushd $MR_ROOT/Libs/External/xgboost
	if [ $OVERRIDE_ALL -gt 0 ]; then
		make clean
		cp $MR_ROOT/Libs/External/xgboost/CMakeLists.txt.orig $MR_ROOT/Libs/External/xgboost/CMakeLists.txt
		cmake .

		#The clean doesn't delete .a files
		find . -name '*.a' -exec rm {} \;
	fi
	make -j20
	cp ./lib/libxgboost.a ${STATIC_LIBS_TARGET}
	cp ./rabit/lib/librabit.a ${STATIC_LIBS_TARGET}
	cp ./dmlc-core/libdmlc.a ${STATIC_LIBS_TARGET}
	popd
fi

if [ "$BUILD_LIGHTGBM" = true ] ; then
	mkdir -p ${BLDDIR}/LightGBM/
	cp ${DIR}/Resources/LightGBM.ubuntu/LightGBM_CMakeLists_file ${BLDDIR}/LightGBM/CMakeLists.txt
	pushd ${BLDDIR}/LightGBM/
	cmake .
	make -j20
	popd 
fi

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
