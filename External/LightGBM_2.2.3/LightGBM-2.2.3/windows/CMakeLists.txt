set(H_FILES 
)
			
set(SRC_FILES
)

# You may create a file named CMakeCustomize.txt to override CMake commands or change flags for a specific project
include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeCustomize.txt OPTIONAL)
if(${CMAKELISTS_OVERRIDE})
  return()
endif()

add_library(windows STATIC ${H_FILES} ${SRC_FILES})


#include("$ENV{MR_ROOT}/Projects/Resources/CMakeUtils/CMakePrintProperties.txt")
