cmake_minimum_required(VERSION 3.5)
project(AlgoMarker VERSION 1.1.0 LANGUAGES C CXX)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/Linux/${CMAKE_BUILD_TYPE}")

#Release Flag settings - optimizations with no debug symbols (faster compile time)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DGIT_HEAD_VERSION='\"'$(GIT_HEAD_VERSION)'\"' -DSO_COMPILATION -Wall -Wno-unused-function -Wno-unknown-pragmas -Wno-comment -Wno-sign-compare -Wno-reorder -Wunused-variable -Wformat -Wparentheses -Wmissing-braces -Wreturn-type -Werror=return-type --std=c++17 -Wno-write-strings -Wuninitialized -fopenmp -msse2 -msse3 -msse4 -march=native --param inline-unit-growth=1000000 -ldl -Ofast -fomit-frame-pointer -fPIC -DHAVE_WORKING_ISFINITE")

#Debug symbols with performance information for tracking program execution times (but still in relase with optimizations, some debug info won't be available)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGIT_HEAD_VERSION='\"'$(GIT_HEAD_VERSION)'\"' -Wall -Wno-unused-function -Wno-unknown-pragmas -Wno-comment -Wno-sign-compare -Wno-reorder -Wunused-variable -Wformat -Wparentheses -Wmissing-braces -Wreturn-type -Werror=return-type --std=c++17 -Wno-write-strings -Wuninitialized -Wreturn-type -fopenmp -msse2 -msse3 -msse4 -march=native --param inline-unit-growth=1000000 -Ofast -pg -g -ldl")

#Debug symbols (but still in relase with optimizations, some debug info won't be available)
#set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} -DGIT_HEAD_VERSION='\"'$(GIT_HEAD_VERSION)'\"' -Wall -Wno-unused-function -Wno-unknown-pragmas -Wno-comment -Wno-sign-compare -Wno-reorder -Wunused-variable -Wformat -Wparentheses -Wmissing-braces -Wreturn-type -Werror=return-type --std=c++17 -Wno-write-strings -Wuninitialized -Wreturn-type -fopenmp -msse2 -msse3 -msse4 -march=native --param inline-unit-growth=1000000 -Ofast -ldl -g")

#for platform independent:
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DGIT_HEAD_VERSION='\"'$(GIT_HEAD_VERSION)'\"' -Wall -Wno-unused-function -Wno-unknown-pragmas -Wno-comment -Wno-sign-compare -Wno-reorder -Wunused-variable -Wformat -Wparentheses -Wmissing-braces -Wreturn-type -Werror=return-type --std=c++17 -Wno-write-strings -Wuninitialized -fopenmp -march=x86-64 -msse2 -msse3 -msse4 --param inline-unit-growth=1000000 -ldl -Ofast -fomit-frame-pointer -fPIC -DHAVE_WORKING_ISFINITE")
#set(CMAKE_EXE_LINKER_FLAGS "-Wl,-Bstatic,--whole-archive -L$ENV{MR_ROOT}/Libs/Internal/AlgoMarker/LinuxSharedLib/Build/lib/ -lboost_regex -lboost_program_options -lboost_system -lboost_filesystem -lxgboost -lrabit -ldmlc -l_lightgbm -Wl,-no-whole-archive -Wl,-Bdynamic -licudata -licui18n -licuuc" CACHE INTERNAL "" FORCE)

# Build with older GlibC
include_directories(SYSTEM /nas1/Work/Libs/Glibc_ubuntu_22.04/include)
include_directories(SYSTEM /nas1/Work/Libs/Glibc_ubuntu_22.04/include/x86_64-linux-gnu)
add_link_options( -L/nas1/Work/Libs/Glibc_ubuntu_22.04)
# End Old glibc

# You may create a file named CMakeCustomize.txt to override CMake commands or change flags for a specific project
#include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeCustomize.txt OPTIONAL)
#if(${CMAKELISTS_OVERRIDE})
#  return()
#endif()

include_directories($ENV{MR_ROOT}/Libs)
include_directories(SYSTEM /nas1/Work/Libs/Boost/boost_1_85_0)
include_directories($ENV{MR_ROOT}/Libs/External)
include_directories($ENV{MR_ROOT}/Libs/Internal)
include_directories(${PROJECT_SOURCE_DIR})
include_directories($ENV{MR_ROOT}/Libs/External/xgboost/src)
include_directories($ENV{MR_ROOT}/Libs/External/xgboost/include)
include_directories($ENV{MR_ROOT}/Libs/External/xgboost/rabit/include)
include_directories($ENV{MR_ROOT}/Libs/External/xgboost/dmlc-core/include)
include_directories($ENV{MR_ROOT}/Libs/External/LightGBM_2.2.3/LightGBM-2.2.3/include)


add_subdirectory($ENV{MR_ROOT}/Libs/Internal/InfraMed/InfraMed InfraMed)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/Logger/Logger Logger)
add_subdirectory($ENV{MR_ROOT}/Libs/External/Mars/Mars Mars)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/MedAlgo/MedAlgo MedAlgo)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/MedEmbed/MedEmbed MedEmbed)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/MedIO/MedIO MedIO)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/MedMat/MedMat MedMat)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/MedPlotly/MedPlotly MedPlotly)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/MedProcessTools/MedProcessTools MedProcessTools)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/MedSparseMat/MedSparseMat MedSparseMat)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/MedSplit/MedSplit MedSplit)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/MedStat/MedStat MedStat)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/MedTime/MedTime MedTime)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/MedUtils/MedUtils MedUtils)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/QRF/QRF QRF)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/SerializableObject/SerializableObject SerializableObject)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/TQRF/TQRF TQRF)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/micNet/micNet micNet)

add_subdirectory($ENV{MR_ROOT}/Libs/Internal/AlgoMarker/DllAPITester DllAPITester)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/AlgoMarker/DynAMWrapper DynAMWrapper)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/AlgoMarker/AMApiTester AMApiTester)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/AlgoMarker/SignalsSumAM SignalsSumAM)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/AlgoMarker/CommonTestingTools CommonTestingTools)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/AlgoMarker/ApplyTool ApplyTool)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/AlgoMarker/StressTest StressTest)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/AlgoMarker/MemroyLeak_Test MemroyLeak_Test)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/AlgoMarker/AlpineCompitible AlpineCompitible)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/AlgoMarker/TestLibSimple TestLibSimple)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/AlgoMarker/InternalAPITester InternalAPITester)
add_subdirectory($ENV{MR_ROOT}/Libs/Internal/AlgoMarker/AlgoMarker AlgoMarker)


