cmake_minimum_required(VERSION 3.5.0)
file(GLOB SRC_FILES
     "*.h"
     "*.cpp"
)

set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} -DHAVE_WORKING_ISFINITE -fvisibility=hidden -Wl,--version-script=$ENV{MR_ROOT}/Libs/Internal/AlgoMarker/AlgoMarker/AlgoMarker-export.version")
add_library(dyn_AlgoMarker SHARED ${SRC_FILES})
add_library(AlgoMarker STATIC ${SRC_FILES})

target_link_libraries(dyn_AlgoMarker -Wl,-Bstatic,--whole-archive InfraMed Logger Mars MedAlgo MedEmbed MedIO MedMat MedProcessTools MedSparseMat MedSplit MedStat MedTime MedUtils QRF SerializableObject TQRF micNet -L/nas1/Work/Libs/Boost/boost_1_85_0/installation/lib/static -lboost_regex -lboost_program_options -lboost_system -lboost_filesystem -lboost_atomic -L/nas1/Work/SharedLibs/linux/ubuntu/static_libs/Release_algomarker -lxgboost -ldmlc -l_lightgbm  -Wl,-no-whole-archive -Wl,-Bdynamic)

