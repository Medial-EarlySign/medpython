cmake_minimum_required(VERSION 3.2)
project(_SOLNAME_)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/Linux/${CMAKE_BUILD_TYPE}")

#Release Flag settings - optimizations with no debug symbols (faster compile time)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DGIT_HEAD_VERSION='\"'$(GIT_HEAD_VERSION)'\"' -Wall -Wno-unused-function -Wno-unknown-pragmas -Wno-comment -Wno-sign-compare -Wno-reorder -Wunused-variable -Wformat -Wparentheses -Wmissing-braces -Wreturn-type -Werror=return-type --std=c++17 -Wno-write-strings -Wuninitialized -fopenmp -msse2 -msse3 -msse4 -march=native --param inline-unit-growth=1000000 -ldl -Ofast -fomit-frame-pointer")

#Debug symbols with performance information for tracking program execution times (but still in relase with optimizations, some debug info won't be available)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGIT_HEAD_VERSION='\"'$(GIT_HEAD_VERSION)'\"' -Wall -Wno-unused-function -Wno-unknown-pragmas -Wno-comment -Wno-sign-compare -Wno-reorder -Wunused-variable -Wformat -Wparentheses -Wmissing-braces -Wreturn-type -Werror=return-type --std=c++17 -Wno-write-strings -Wuninitialized -Wreturn-type -fopenmp -msse2 -msse3 -msse4 -march=native --param inline-unit-growth=1000000 -Ofast -pg -g -ldl")

#Debug symbols (but still in relase with optimizations, some debug info won't be available)
#set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} -DGIT_HEAD_VERSION='\"'$(GIT_HEAD_VERSION)'\"' -Wall -Wno-unused-function -Wno-unknown-pragmas -Wno-comment -Wno-sign-compare -Wno-reorder -Wunused-variable -Wformat -Wparentheses -Wmissing-braces -Wreturn-type -Werror=return-type --std=c++17 -Wno-write-strings -Wuninitialized -Wreturn-type -fopenmp -msse2 -msse3 -msse4 -march=native --param inline-unit-growth=1000000 -Ofast -ldl -g")

#for platform independent:
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DGIT_HEAD_VERSION='\"'$(GIT_HEAD_VERSION)'\"' -Wall -Wno-unused-function -Wno-unknown-pragmas -Wno-comment -Wno-sign-compare -Wno-reorder -Wunused-variable -Wformat -Wparentheses -Wmissing-braces -Wreturn-type -Werror=return-type --std=c++17 -Wno-write-strings -Wuninitialized -fopenmp -march=x86-64 -msse2 -msse3 -msse4 --param inline-unit-growth=1000000 -ldl -Ofast -fomit-frame-pointer -fPIC -DHAVE_WORKING_ISFINITE")
#set(CMAKE_EXE_LINKER_FLAGS "-Wl,-Bstatic,--whole-archive -L$ENV{MR_ROOT}/Libs/Internal/AlgoMarker/LinuxSharedLib/Build/lib/ -lboost_regex -lboost_program_options -lboost_system -lboost_filesystem -lxgboost -lrabit -ldmlc -l_lightgbm -Wl,-no-whole-archive -Wl,-Bdynamic -licudata -licui18n -licuuc" CACHE INTERNAL "" FORCE)

include("$ENV{MR_ROOT}/Projects/Resources/CMakeUtils/CMakeDropNDEBUG.txt")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -march=native -Og")

# You may create a file named CMakeCustomize.txt to override CMake commands or change flags for a specific project
include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeCustomize.txt OPTIONAL)
if(${CMAKELISTS_OVERRIDE})
  return()
endif()

include_directories($ENV{MR_ROOT}/Libs)
include_directories(SYSTEM /server/Work/Libs/Boost/boost_1_67_0-fPIC.ubuntu)
include_directories($ENV{MR_ROOT}/Libs/External)
include_directories($ENV{MR_ROOT}/Libs/Internal)
include_directories($ENV{MR_ROOT}/Libs/InternalOrig)
include_directories(${PROJECT_SOURCE_DIR})
include_directories($ENV{MR_ROOT}/Libs/External/xgboost/src)
include_directories($ENV{MR_ROOT}/Libs/External/xgboost/include)
include_directories($ENV{MR_ROOT}/Libs/External/xgboost/rabit/include)
include_directories($ENV{MR_ROOT}/Libs/External/xgboost/dmlc-core/include)
include_directories($ENV{MR_ROOT}/Libs/External/LightGBM_2.2.3/LightGBM-2.2.3/include)

add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/AlgoMarker/AlgoMarker AlgoMarker)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/AlgoMarker/InternalAPITester InternalAPITester)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/InfraMed/InfraMed InfraMed)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/Logger/Logger Logger)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/External/Mars/Mars Mars)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/MedAlgo/MedAlgo MedAlgo)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/MedEmbed/MedEmbed MedEmbed)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/MedIO/MedIO MedIO)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/MedMat/MedMat MedMat)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/MedPlotly/MedPlotly MedPlotly)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/MedProcessTools/MedProcessTools MedProcessTools)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/MedSparseMat/MedSparseMat MedSparseMat)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/MedSplit/MedSplit MedSplit)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/MedStat/MedStat MedStat)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/MedTime/MedTime MedTime)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/MedUtils/MedUtils MedUtils)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/QRF/QRF QRF)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/SerializableObject/SerializableObject SerializableObject)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/TQRF/TQRF TQRF)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/micNet/micNet micNet)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/AlgoMarker/DllAPITester DllAPITester)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/AlgoMarker/DynAMWrapper DynAMWrapper)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/AlgoMarker/AMApiTester AMApiTester)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/AlgoMarker/SignalsSumAM SignalsSumAM)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/AlgoMarker/CommonTestingTools CommonTestingTools)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/AlgoMarker/ApplyTool ApplyTool)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/AlgoMarker/StressTest StressTest)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/AlgoMarker/MemroyLeak_Test MemroyLeak_Test)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/AlgoMarker/AlpineCompitible AlpineCompitible)
add_subdirectory(/nas1/UsersData/alon/MR/Libs/Internal/AlgoMarker/TestLibSimple TestLibSimple)


