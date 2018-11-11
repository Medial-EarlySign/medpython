#!/bin/bash

RELEASE_PATH=${MR_ROOT}/Libs/Internal/MedPyExport/generate_binding/CMakeBuild/Linux/Release/MedPython
NEW_RELEASE_PATH=${MR_ROOT}/Libs/Internal/MedPyExport/generate_binding/Release

#medial python 3
export PYTHON_INCLUDE_DIR=/opt/medial/dist/usr/include/python3.6m
export PYTHON_LIBRARY=/opt/medial/dist/usr/lib/libpython3.so
if [ -f $PYTHON_LIBRARY ]; then
${MR_ROOT}/Libs/Internal/MedPyExport/generate_binding/make.sh
fi

#anaconda 2
export PYTHON_INCLUDE_DIR=/home/python/anaconda2/include/python2.7
export PYTHON_LIBRARY=/home/python/anaconda2/lib/libpython2.7.so
if [ -f $PYTHON_LIBRARY ]; then
${MR_ROOT}/Libs/Internal/MedPyExport/generate_binding/make.sh
fi

#medial python 2
export PYTHON_INCLUDE_DIR=/opt/medial/python27/usr/include/python2.7
export PYTHON_LIBRARY=/opt/medial/python27/usr/lib/libpython2.7.so
if [ -f $PYTHON_LIBRARY ]; then
${MR_ROOT}/Libs/Internal/MedPyExport/generate_binding/make.sh
fi

#centos RH
export PYTHON_INCLUDE_DIR=/usr/include/python2.7
export PYTHON_LIBRARY=/usr/lib64/libpython2.7.so
if [ -f $PYTHON_LIBRARY ]; then
${MR_ROOT}/Libs/Internal/MedPyExport/generate_binding/make.sh
fi

