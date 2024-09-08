cmake_minimum_required(VERSION 3.5.0)
file(GLOB SRC_FILES
     "*.h"
     "*.cpp"
)
add_library(AlpineCompitible SHARED ${SRC_FILES})