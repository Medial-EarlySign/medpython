#!/bin/bash
BLDDIR=${0%/*}
pushd ${BLDDIR}

version_txt=`get_git_status_text.py | sed 's|"||g' | awk -F"\t" '{print "\"" $1 "\""}'`
echo -e "Git version info:\n${version_txt}"
touch ${MR_ROOT}/Libs/Internal/MedUtils/MedUtils/MedGitVersion.h

pushd ${BLDDIR}/build/Release
make -j20 -e GIT_HEAD_VERSION="$version_txt"
popd 
cp ${BLDDIR}/build/Release/AlgoMarker/libdyn_AlgoMarker.so ${BLDDIR}/Linux/Release/
strip ${BLDDIR}/Linux/Release/libdyn_AlgoMarker.so
nm -CD ${BLDDIR}/Linux/Release/libdyn_AlgoMarker.so | grep " T " > ${BLDDIR}/Linux/Release/AlgoMarker.public_symbols
ldd ${BLDDIR}/Linux/Release/libdyn_AlgoMarker.so > ${BLDDIR}/Linux/Release/AlgoMarker.ldd
FINAL_PT=$(realpath ${BLDDIR}/Linux/Release/libdyn_AlgoMarker.so)
echo "AlgoMarker Shared Object file written to ${FINAL_PT}"
popd