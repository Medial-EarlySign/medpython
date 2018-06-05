project(LightGBM)
cmake_minimum_required(VERSION 2.8)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/Linux/${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} --std=c++11 -Wno-write-strings -Wuninitialized -fopenmp -msse2 -msse3 -msse4 -march=native --param inline-unit-growth=1000000 -fomit-frame-pointer -Ofast")
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} --std=c++11 -Wno-write-strings -Wuninitialized -fopenmp -msse2 -msse3 -msse4 -march=native --param inline-unit-growth=1000000 -Ofast -pg")
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} --std=c++11 -Wno-write-strings -Wuninitialized -fopenmp -msse2 -msse3 -msse4 -march=native -O0 -g")
include("$ENV{MR_ROOT}/Projects/Resources/CMakeUtils/CMakeDropNDEBUG.txt")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -march=native -Og")

include_directories($ENV{MR_ROOT}/Libs)
include_directories($ENV{MR_ROOT}/Libs/External)
include_directories($ENV{MR_ROOT}/Libs/Internal)
include_directories($ENV{MR_ROOT}/Libs/InternalOrig)
include_directories($ENV{R_INCL_DIR})
include_directories($ENV{LIBXL_INCL_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories($ENV{MR_ROOT}/Libs/External/xgboost/src)
include_directories($ENV{MR_ROOT}/Libs/External/xgboost/include)
include_directories($ENV{MR_ROOT}/Libs/External/xgboost/rabit/include)
include_directories($ENV{MR_ROOT}/Libs/External/xgboost/dmlc-core/include)
include_directories($ENV{MR_ROOT}/Libs/External/LightGBM/LightGBM/include)

add_subdirectory(windows)


#include("$ENV{MR_ROOT}/Projects/Resources/CMakeUtils/CMakePrintProperties.txt")
