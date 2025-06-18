cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0086 NEW)
cmake_policy(SET CMP0078 NEW)

#execute_process(COMMAND swig -python medpython.i
#  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

#This always finds the system's python installtion
#FIND_PACKAGE(PythonLibs)
message("*********************************************")
if("${Python3_INCLUDE_DIRS}" STREQUAL "")
    message(FATAL_ERROR "Could not find PYTHON_INCLUDE_DIR")
endif()
if("${Python3_LIBRARIES}" STREQUAL "")
    message(FATAL_ERROR "Could not find Python3_LIBRARIES")
endif()
message("*********************************************")

ADD_CUSTOM_TARGET(apply_directives
                  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/make_apply.py ${CMAKE_CURRENT_SOURCE_DIR}/../../MedPyExport/ ${CMAKE_CURRENT_SOURCE_DIR}/apply_directives.i )


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${Python3_INCLUDE_DIRS})

SET(CMAKE_SWIG_FLAGS "")
SET_SOURCE_FILES_PROPERTIES(medpython.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(medpython.i PROPERTIES SWIG_FLAGS "-includeall")
message("\n*********************************************\n${CMAKE_CXX_FLAGS}\n*********************************************\n")
file(GLOB CppFiles ${CMAKE_CURRENT_SOURCE_DIR}/../../MedPyExport/*.cpp)
file(GLOB HeaderFiles ${CMAKE_CURRENT_SOURCE_DIR}/../../MedPyExport/*.h)


SET_SOURCE_FILES_PROPERTIES(SOURCE ${CppFiles} PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -lpython3 -lpthread -ldl -lutil -lm")
set(SWIG_MODULE_medpython_EXTRA_DEPS MedPython.h MedPython.h ${HeaderFiles} apply_directives pythoncode.i)

set_property(SOURCE medpython.i PROPERTY SWIG_MODULE_NAME medpython)

swig_add_library(medpython LANGUAGE python SOURCES medpython.i ${CppFiles})


SWIG_LINK_LIBRARIES(medpython ${Python3_LIBRARIES} -Wl,-Bstatic,--whole-archive -Wl,--start-group ${MEDIAL_INTERNAL_LIBS} Mars ${BOOST_STATIC_LIB_PATHS} xgboost _lightgbm_static -Wl,--end-group -Wl,-no-whole-archive -Wl,-Bdynamic)