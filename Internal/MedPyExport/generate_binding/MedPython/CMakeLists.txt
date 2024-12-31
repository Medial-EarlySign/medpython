cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0086 NEW)
cmake_policy(SET CMP0078 NEW)

#execute_process(COMMAND swig -python medpython.i
#  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

#This always finds the system's python installtion
#FIND_PACKAGE(PythonLibs)
#message("cmake found PYTHON_INCLUDE_PATH: ${PYTHON_INCLUDE_PATH}")
#message("cmake found PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")

message("*********************************************")
if("$ENV{PYTHON_INCLUDE_DIR}" STREQUAL "")
    message(FATAL_ERROR "Could not find PYTHON_INCLUDE_DIR environment variable")
endif()
if("$ENV{PYTHON_LIBRARY}" STREQUAL "")
    message(FATAL_ERROR "Could not find PYTHON_LIBRARY environment variable")
endif()
message("PYTHON_INCLUDE_DIR=$ENV{PYTHON_INCLUDE_DIR}")
message("PYTHON_LIBRARY=$ENV{PYTHON_LIBRARY}")
message("*********************************************")

ADD_CUSTOM_TARGET(apply_directives
                  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/make_apply.py ${CMAKE_CURRENT_SOURCE_DIR}/../../MedPyExport/ ${CMAKE_CURRENT_SOURCE_DIR}/apply_directives.i )


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES($ENV{PYTHON_INCLUDE_DIR})

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


SWIG_LINK_LIBRARIES(medpython $ENV{PYTHON_LIBRARY} -Wl,-Bstatic,--whole-archive MedMat MedSplit InfraMed Logger Mars MedAlgo MedEmbed MedFeat MedIO MedProcessTools MedSparseMat MedStat MedTime MedUtils QRF SerializableObject TQRF micNet CommonLib  /server/Work/Libs/Boost/boost_1_85_0/installation/lib/static/libboost_regex.a /server/Work/Libs/Boost/boost_1_85_0/installation/lib/static/libboost_program_options.a /server/Work/Libs/Boost/boost_1_85_0/installation/lib/static/libboost_filesystem.a /server/Work/Libs/Boost/boost_1_85_0/installation/lib/static/libboost_system.a /server/Work/Libs/Boost/boost_1_85_0/installation/lib/static/libboost_atomic.a /nas1/Work/SharedLibs/linux/ubuntu/static_libs/Release_new/libxgboost.a /nas1/Work/SharedLibs/linux/ubuntu/static_libs/Release_new/lib_lightgbm.a  /nas1/Work/SharedLibs/linux/ubuntu/static_libs/Release_new/libdmlc.a -Wl,-no-whole-archive -Wl,-Bdynamic)

#
#  TODO: Why doesn't it compile with: 
#  $ENV{BOOST_PO_LIB_${CMAKE_BUILD_TYPE}} $ENV{BOOST_FS_LIB_${CMAKE_BUILD_TYPE}} $ENV{BOOST_RE_LIB_${CMAKE_BUILD_TYPE}} $ENV{BOOST_SYS_LIB_${CMAKE_BUILD_TYPE}}
#
#


#set(H_FILES 
#)
##
##set(SRC_FILES
##    MedPython.cpp
##    MedPython_wrap.c
##)
##
##if(NOT DEFINED $ENV{PYTHON_INCLUDE_DIRS})
##    message("Please define PYTHON_INCLUDE_DIRS or use flag -DPYTHON_INCLUDE_DIRS")
##endif()
##
##include_directories($ENV{PYTHON_INCLUDE_DIRS})
##
##
##add_library(MedPython SHARED ${SRC_FILES})
##
##set_target_properties(MedPython PROPERTIES PREFIX "")
##set_target_properties(MedPython PROPERTIES LIBRARY_OUTPUT_NAME _medpython)

#include("$ENV{MR_ROOT}/Projects/Resources/CMakeUtils/CMakePrintProperties.txt")
