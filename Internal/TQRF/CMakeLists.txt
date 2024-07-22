project(TQRF)
cmake_minimum_required(VERSION 2.8)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/Linux/${CMAKE_BUILD_TYPE}")

#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} --std=c++11 -Wno-write-strings -Wuninitialized -fopenmp -msse2 -msse3 -msse4 -march=native --param inline-unit-growth=1000000 -Ofast -g")
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} --std=c++11 -Wno-write-strings -Wuninitialized -fopenmp -msse2 -msse3 -msse4 -march=native --param inline-unit-growth=1000000 -Ofast -pg")
set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} --std=c++11 -Wno-write-strings -Wuninitialized -fopenmp -msse2 -msse3 -msse4 -march=native --param inline-unit-growth=1000000 -Ofast -fomit-frame-pointer")

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

add_subdirectory(TQRF)
add_subdirectory(TQRF_Tester)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/InfraMed/InfraMed InfraMed)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/Logger/Logger Logger)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/External/Mars/Mars Mars)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/MedAlgo/MedAlgo MedAlgo)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/MedFeat/MedFeat MedFeat)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/MedProcessTools/MedProcessTools MedProcessTools)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/MedStat/MedStat MedStat)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/MedTime/MedTime MedTime)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/MedUtils/MedUtils MedUtils)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/QRF/QRF QRF)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/External/gbm/gbm gbm)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/micNet/micNet micNet)


#include("$ENV{MR_ROOT}/Projects/Resources/CMakeUtils/CMakePrintProperties.txt")
