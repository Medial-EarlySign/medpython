cmake_minimum_required(VERSION 3.14)
project(AlgoMarker_Tools VERSION 1.1.0 LANGUAGES CXX)
cmake_policy(SET CMP0074 NEW)
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.27")
    cmake_policy(SET CMP0144 NEW)
endif()

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/Linux/${CMAKE_BUILD_TYPE}")

if(UNIX AND NOT APPLE)
    include(CheckIncludeFileCXX)
    
    # Check if execinfo.h exists (Standard Linux/glibc has it, Alpine/Musl does not)
    check_include_file_cxx("execinfo.h" HAVE_EXECINFO)

    if(NOT HAVE_EXECINFO)
        message(STATUS "Musl libc (Alpine) detected. Applying compatibility flags.")
        
        # Add the flags specifically for this environment
        add_compile_definitions(
            DMLC_LOG_STACK_TRACE=0
            DMLC_CMAKE_LITTLE_ENDIAN=1
            mmap64=mmap
        )
    endif()
endif()

if(DEFINED ENV{BOOST_ROOT})
    message(STATUS "BOOST_ROOT is set to: $ENV{BOOST_ROOT}")
    set(BOOST_ROOT "$ENV{BOOST_ROOT}")
else()
    message(STATUS "BOOST_ROOT environment variable is not set. Please set it to the Boost installation path.")
    # set(BOOST_ROOT "$ENV{HOME}/Documents/MES/Boost") #If you have compiled boost by yourself with -fPIC flag, please provide the path to folder with "/lib" with shared and static lib and "/include" for headers.
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -fopenmp -Wall -Wno-write-strings -Wno-sign-compare -Wno-reorder -Wno-unknown-pragmas -Wno-unused-function -ldl  --param inline-unit-growth=1000000 -msse4 -DHAVE_WORKING_ISFINITE -DSO_COMPILATION -DGIT_HEAD_VERSION='\"'$(GIT_HEAD_VERSION)'\"'")

set(BOOST_LIBS
    regex
    filesystem
    system
    program_options
) # for certain version, you might need to add "atomic"

find_package(Boost REQUIRED)
if(Boost_VERSION VERSION_GREATER_EQUAL "1.69") # 1.69 is represented as 106900 in CMake
    # Remove "system" for newer versions (>= 1.69)
    # The list(REMOVE_ITEM) command is used to remove a specific item from a list variable.
    list(REMOVE_ITEM BOOST_LIBS "system")
    message(STATUS "Boost version is ${Boost_VERSION} (>= 1.69), removing 'system' from BOOST_LIBS.")
else()
    message(STATUS "Boost version is ${Boost_VERSION} (< 1.69), keeping 'system' in BOOST_LIBS.")
endif()

set (MEDIAL_INTERNAL_LIBS InfraMed Logger MedAlgo MedEmbed MedIO MedMat MedPlotly MedProcessTools MedSparseMat MedSplit MedStat MedTime MedUtils QRF SerializableObject TQRF micNet
) #CommonLib if needed

# Dynamic linking
function(add_linking_flags TARGET_NAME)
    set(BUILD_STATIC_LIB OFF CACHE BOOL "Build XGBoost as a dynamic library" FORCE)
    if (BOOST_ROOT)
        set(Boost_LIBRARIES "")
        foreach(lib IN LISTS BOOST_LIBS)
            list(APPEND Boost_LIBRARIES ${BOOST_ROOT}/lib/libboost_${lib}.so)
        endforeach()
    endif()
    target_link_libraries(${TARGET_NAME} ${ADDITIONAL_LINK_FLAGS} -Wl,--start-group ${MEDIAL_INTERNAL_LIBS} DynAMWrapper AlgoMarker CommonTestingTools xgboost _lightgbm  -Wl,--end-group  ${Boost_LIBRARIES} )
endfunction()

#Static linking
function(add_linking_flags_static TARGET_NAME)
    set(BOOST_STATIC_LIB_PATHS "")
    set(Boost_USE_STATIC_LIBS ON)
    foreach(lib IN LISTS BOOST_LIBS)
	    string(TOUPPER "${lib}" LIB_VAR)

        if (NOT BOOST_ROOT)
            find_library(BOOST_${LIB_VAR}_STATIC
                NAMES "libboost_${lib}.a"
            )
        else()
            set(BOOST_${LIB_VAR}_STATIC ${BOOST_ROOT}/lib/libboost_${lib}.a)
        endif()

        if(NOT BOOST_${LIB_VAR}_STATIC)
            message(FATAL_ERROR "Could not find static Boost library: boost_${lib}.a")
        endif()

        list(APPEND BOOST_STATIC_LIB_PATHS ${BOOST_${LIB_VAR}_STATIC})
    endforeach()

    set(BUILD_STATIC_LIB ON CACHE BOOL "Build XGBoost as a static library" FORCE)
    target_link_libraries(${TARGET_NAME} ${ADDITIONAL_LINK_FLAGS} -Bstatic,--whole-archive -Wl,--start-group ${MEDIAL_INTERNAL_LIBS} -Wl,--end-group xgboost _lightgbm_static  ${BOOST_STATIC_LIB_PATHS} -Wl,-no-whole-archive -Wl,-Bdynamic)

endfunction()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../External)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal)
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../External/xgboost/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../External/xgboost/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../External/xgboost/rabit/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../External/xgboost/dmlc-core/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../External/LightGBM_2.2.3/LightGBM-2.2.3/include)

# Find Boost
find_package(OpenMP REQUIRED)

find_package(Boost REQUIRED COMPONENTS ${BOOST_LIBS})  # Example components

if (BOOST_ROOT)
    set(Boost_INCLUDE_DIRS "${BOOST_ROOT}/include")
endif()
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
# Add Internal libraries
foreach(mes_lib IN LISTS MEDIAL_INTERNAL_LIBS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/${mes_lib}/${mes_lib} ${mes_lib})
endforeach()

# External Libraries
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../External/xgboost xgboost)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../External/LightGBM_2.2.3/LightGBM-2.2.3 _lightgbm)


# Or add each project:
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/DllAPITester DllAPITester)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/DynAMWrapper DynAMWrapper)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/CommonTestingTools CommonTestingTools)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/AMApiTester AMApiTester)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/StressTest StressTest)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/MemroyLeak_Test MemroyLeak_Test)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/TestLibSimple TestLibSimple)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/InternalAPITester InternalAPITester)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/AlgoMarker AlgoMarker)

# Inside each project:
#file(GLOB SRC_FILES
#     "*.h"
#     "*.cpp"
#)
#add_executable(Test ${SRC_FILES})
#add_linking_flags_dynamic(Test)
#add_linking_flags_static(Test)
