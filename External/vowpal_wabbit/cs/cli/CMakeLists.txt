set(H_FILES 
	clr_io.h
	clr_io_memory.h
	resource.h
	spanning_tree_clr.h
	vld_clr.h
	vowpalwabbit.h
	vw_arguments.h
	vw_base.h
	vw_builder.h
	vw_cbutil.h
	vw_clr.h
	vw_example.h
	vw_interface.h
	vw_labelcomparator.h
	vw_label.h
	vw_model.h
	vw_prediction.h
	vw_settings.h
)
			
set(SRC_FILES
	AssemblyInfo.cpp
	clr_io.cpp
	clr_io_memory.cpp
	spanning_tree_clr.cpp
	vld_clr.cpp
	vowpalwabbit.cpp
	vw_base.cpp
	vw_builder.cpp
	vw_cbutil.cpp
	vw_example.cpp
	vw_exception.cpp
	vw_model.cpp
	vw_prediction.cpp
)

# You may create a file named CMakeCustomize.txt to override CMake commands or change flags for a specific project
include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeCustomize.txt OPTIONAL)
if(${CMAKELISTS_OVERRIDE})
  return()
endif()

add_library(cli STATIC ${H_FILES} ${SRC_FILES})
add_library(dyn_cli SHARED ${H_FILES} ${SRC_FILES})


#include("$ENV{MR_ROOT}/Projects/Resources/CMakeUtils/CMakePrintProperties.txt")
