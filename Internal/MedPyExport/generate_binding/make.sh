cp SWIG.CMakeLists.txt CMakeLists.txt
cp MedPython/SWIG.CMakeLists.txt MedPython/CMakeLists.txt
mkdir -p $MR_ROOT/Libs/Internal/MedPyExport/generate_binding/CMakeBuild/Linux/Release
pushd $MR_ROOT/Libs/Internal/MedPyExport/generate_binding/CMakeBuild/Linux/Release 
cmake ../../../
make -j 8;
popd
