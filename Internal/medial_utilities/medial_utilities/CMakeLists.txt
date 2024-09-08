set(H_FILES 
	globalRNG.h
	medial_utilities.h
)
			
set(SRC_FILES
	clustering.cpp
	io.cpp
	numerical.cpp
	utilities.cpp
	vector_utils.cpp
	vector_utils_dbl.cpp
)

# You may create a file named CMakeCustomize.txt to override CMake commands or change flags for a specific project
include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeCustomize.txt OPTIONAL)
if(${CMAKELISTS_OVERRIDE})
  return()
endif()

add_library(medial_utilities STATIC ${H_FILES} ${SRC_FILES})


#include("$ENV{MR_ROOT}/Projects/Resources/CMakeUtils/CMakePrintProperties.txt")
