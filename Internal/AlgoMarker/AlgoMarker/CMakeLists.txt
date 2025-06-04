cmake_minimum_required(VERSION 3.5.0)
file(GLOB SRC_FILES
     "*.h"
     "*.cpp"
)

set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} -DHAVE_WORKING_ISFINITE -fvisibility=hidden -Wl,--version-script=${CMAKE_CURRENT_SOURCE_DIR}/AlgoMarker-export.version")
add_library(dyn_AlgoMarker SHARED ${SRC_FILES})
add_library(AlgoMarker STATIC ${SRC_FILES})

add_linking_flags_static(dyn_AlgoMarker)

