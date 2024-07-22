project(MedPyExport)
cmake_minimum_required(VERSION 2.8)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/Linux/${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DGIT_HEAD_VERSION='\"'$(GIT_HEAD_VERSION)'\"' -Wall -Wno-unused-function -Wno-unknown-pragmas -Wno-comment -Wno-sign-compare -Wno-reorder -Wunused-variable -Wformat -Wparentheses -Wmissing-braces -Wreturn-type -Werror=return-type --std=c++11 -Wno-write-strings -Wuninitialized -fopenmp -msse2 -msse3 -msse4 -march=native --param inline-unit-growth=1000000 -ldl -Ofast -g")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGIT_HEAD_VERSION='\"'$(GIT_HEAD_VERSION)'\"' -Wall -Wno-unused-function -Wno-unknown-pragmas -Wno-comment -Wno-sign-compare -Wno-reorder -Wunused-variable -Wformat -Wparentheses -Wmissing-braces -Wreturn-type -Werror=return-type --std=c++11 -Wno-write-strings -Wuninitialized -Wreturn-type -fopenmp -msse2 -msse3 -msse4 -march=native --param inline-unit-growth=1000000 -Ofast -pg -g -ldl")
#set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} -DGIT_HEAD_VERSION='\"'$(GIT_HEAD_VERSION)'\"' -Wall -Wno-unused-function -Wno-unknown-pragmas -Wno-comment -Wno-sign-compare -Wno-reorder -Wunused-variable -Wformat -Wparentheses -Wmissing-braces -Wreturn-type -Werror=return-type --std=c++11 -Wno-write-strings -Wuninitialized -Wreturn-type -fopenmp -msse2 -msse3 -msse4 -march=native --param inline-unit-growth=1000000 -Ofast -fomit-frame-pointer -ldl")

include("$ENV{MR_ROOT}/Projects/Resources/CMakeUtils/CMakeDropNDEBUG.txt")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -march=native -Og")

# You may create a file named CMakeCustomize.txt to override CMake commands or change flags for a specific project
include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeCustomize.txt OPTIONAL)
if(${CMAKELISTS_OVERRIDE})
  return()
endif()

include_directories($ENV{MR_ROOT}/Libs)
include_directories($ENV{MR_ROOT}/Libs/External)
include_directories($ENV{MR_ROOT}/Libs/Internal)
include_directories($ENV{MR_ROOT}/Libs/InternalOrig)
include_directories($ENV{R_INCL_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories($ENV{MR_ROOT}/Libs/External/xgboost/src)
include_directories($ENV{MR_ROOT}/Libs/External/xgboost/include)
include_directories($ENV{MR_ROOT}/Libs/External/xgboost/rabit/include)
include_directories($ENV{MR_ROOT}/Libs/External/xgboost/dmlc-core/include)
include_directories($ENV{MR_ROOT}/Libs/External/LightGBM_2.2.3/LightGBM-2.2.3/include)

add_subdirectory(MedPyExport)


#include("$ENV{MR_ROOT}/Projects/Resources/CMakeUtils/CMakePrintProperties.txt")
