cmake_minimum_required(VERSION 3.5.0)
file(GLOB SRC_FILES
     "*.h"
     "*.cpp"
)

add_executable(TestLibSimple ${SRC_FILES})
target_link_libraries(TestLibSimple ${ADDITIONAL_LINK_FLAGS})
