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
    set(Boost_INCLUDE_DIRS "${BOOST_ROOT}/include")
    include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
    set(BUILD_EXEC ON)
else()
    message(STATUS "BOOST_ROOT environment variable is not set. Please set it to the Boost installation path.")
    # set(BOOST_ROOT "$ENV{HOME}/Documents/MES/Boost") #If you have compiled boost by yourself with -fPIC flag, please provide the path to folder with "/lib" with shared and static lib and "/include" for headers.
    find_package(Boost QUIET)
    if (NOT Boost_FOUND)
        message(STATUS "BOOST_ROOT environment variable is not set. Please set it to the Boost installation path, if you want to use a specific boost.")
        message(STATUS "Boost not found. Downloading headers...")
        set(BUILD_EXEC OFF)
        include(FetchContent)
        FetchContent_Declare(
            boost_headers
                URL      https://archives.boost.io/release/1.89.0/source/boost_1_89_0.tar.gz
                SOURCE_SUBDIR pathThatDoesNotExist
                DOWNLOAD_EXTRACT_TIMESTAMP ON
            )
            
            FetchContent_GetProperties(boost_headers)
            if(NOT boost_headers_POPULATED)
                FetchContent_MakeAvailable(boost_headers)
            endif()

            include_directories(${boost_headers_SOURCE_DIR})
    else()
        set(BUILD_EXEC ON)
    endif()
endif()

## Eigen
find_package(Eigen3 QUIET)

if(Eigen3_FOUND)
    message(STATUS "Eigen3 found. Using system headers at: ${EIGEN3_INCLUDE_DIR}")
    include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
else()
    message(STATUS "Eigen3 not found. Downloading headers from source...")
    include(FetchContent)
    FetchContent_Declare(
        eigen_headers
        URL      https://gitlab.com/libeigen/eigen/-/archive/3.4.1/eigen-3.4.1.tar.gz
        SOURCE_SUBDIR pathThatDoesNotExist
        DOWNLOAD_EXTRACT_TIMESTAMP ON
    )
    
    FetchContent_GetProperties(eigen_headers)
    if(NOT eigen_headers_POPULATED)
        FetchContent_MakeAvailable(eigen_headers)
    endif()

    include_directories(SYSTEM ${eigen_headers_SOURCE_DIR})
endif()

# Json:
find_package(nlohmann_json QUIET)
if (nlohmann_json_FOUND)
    message(STATUS "nlohmann_json found. Using system headers at: ${nlohmann_json_INCLUDE_DIR}")
    include_directories(SYSTEM ${nlohmann_json_INCLUDE_DIR})
else ()
    # https://github.com/nlohmann/json/releases/download/v3.11.2/json.hpp include.zip
    include(FetchContent)
    FetchContent_Declare(
        json_headers
        URL      https://github.com/nlohmann/json/releases/download/v3.11.2/include.zip
        SOURCE_SUBDIR pathThatDoesNotExist
        DOWNLOAD_EXTRACT_TIMESTAMP ON
    )
    FetchContent_GetProperties(json_headers)
    if(NOT json_headers_POPULATED)
        FetchContent_MakeAvailable(json_headers)
    endif()
    include_directories(SYSTEM "${json_headers_SOURCE_DIR}/single_include")
endif()

# XGboost
set(XGBOOST_VERSION "2.0.3")

if(XGBOOST_VERSION VERSION_GREATER "2.0.3")
    set(XGBOOST_TARGET_DIR "")
else()
    set(XGBOOST_TARGET_DIR "cpp_src")
endif()

if(XGBOOST_VERSION VERSION_GREATER_EQUAL "3.0.0")
    set(XGBOOST_DOWNLOAD "-src")
else()
    set(XGBOOST_DOWNLOAD "")
endif()

FetchContent_Declare(
    xgboost
    URL      https://github.com/dmlc/xgboost/releases/download/v${XGBOOST_VERSION}/xgboost${XGBOOST_DOWNLOAD}-${XGBOOST_VERSION}.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    SOURCE_SUBDIR ${XGBOOST_TARGET_DIR}
)

set(BUILD_STATIC_LIB ON CACHE BOOL "Build XGBoost as a static library" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static libs" FORCE)
set(XGBOOST_BUILD_TESTS OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(xgboost)



include_directories(SYSTEM ${xgboost_SOURCE_DIR}/${XGBOOST_TARGET_DIR}/include)          
include_directories(SYSTEM ${xgboost_SOURCE_DIR}/${XGBOOST_TARGET_DIR}/dmlc-core/include)

# LightGBM
FetchContent_Declare(
    lightgbm
    URL      https://github.com/lightgbm-org/LightGBM/releases/download/v4.6.0/lightgbm-4.6.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP ON
)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(BUILD_CLI OFF CACHE INTERNAL "")
FetchContent_MakeAvailable(lightgbm)

include_directories(SYSTEM ${lightgbm_SOURCE_DIR}/include)     
include_directories(SYSTEM ${lightgbm_SOURCE_DIR}/external_libs/fast_double_parser/include) 
include_directories(SYSTEM ${lightgbm_SOURCE_DIR}/external_libs/fmt/include) 
## End

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -fopenmp -DMES_LIBRARY=1 -Wall -Wno-write-strings -Wno-sign-compare -Wno-reorder -Wno-unknown-pragmas -Wno-unused-function -ldl  --param inline-unit-growth=1000000 -msse4 -DHAVE_WORKING_ISFINITE -DSO_COMPILATION -DGIT_HEAD_VERSION='\"'$(GIT_HEAD_VERSION)'\"'")

set(BOOST_LIBS
    system
    program_options
) # for certain version, you might need to add "atomic"

find_package(Boost QUIET)
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

    find_package(Boost REQUIRED COMPONENTS ${BOOST_LIBS})  # Example components

    if (BOOST_ROOT)
        set(Boost_LIBRARIES "")
        foreach(lib IN LISTS BOOST_LIBS)
            list(APPEND Boost_LIBRARIES ${BOOST_ROOT}/lib/libboost_${lib}.so)
        endforeach()
    endif()
    target_link_libraries(${TARGET_NAME} ${ADDITIONAL_LINK_FLAGS} -Wl,--start-group ${MEDIAL_INTERNAL_LIBS} DynAMWrapper AlgoMarker CommonTestingTools xgboost _lightgbm  -Wl,--end-group  ${Boost_LIBRARIES} )
endfunction()

#Static linking
set(BUILD_STATIC_LIB ON CACHE BOOL "Build XGBoost as a static library" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static libs" FORCE)

function(add_linking_flags_static TARGET_NAME)
    target_link_libraries(${TARGET_NAME} ${ADDITIONAL_LINK_FLAGS} -Bstatic,--whole-archive -Wl,--start-group ${MEDIAL_INTERNAL_LIBS} -Wl,--end-group xgboost _lightgbm  -Wl,-no-whole-archive -Wl,-Bdynamic)
endfunction()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal)
include_directories(${PROJECT_SOURCE_DIR})

# Find Boost
find_package(OpenMP REQUIRED)



# Add Internal libraries
foreach(mes_lib IN LISTS MEDIAL_INTERNAL_LIBS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/${mes_lib}/${mes_lib} ${mes_lib})
endforeach()

# Or add each project:
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/AlgoMarker AlgoMarker)

if (BUILD_EXEC)
    message(STATUS "BUILD Also Executables for Testing")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/DllAPITester DllAPITester)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/DynAMWrapper DynAMWrapper)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/CommonTestingTools CommonTestingTools)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/AMApiTester AMApiTester)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/StressTest StressTest)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/MemroyLeak_Test MemroyLeak_Test)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/TestLibSimple TestLibSimple)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../Internal/AlgoMarker/InternalAPITester InternalAPITester)
else()
    message(STATUS "SKIP building executables, Boost is missing")
endif()

# Inside each project:
#file(GLOB SRC_FILES
#     "*.h"
#     "*.cpp"
#)
#add_executable(Test ${SRC_FILES})
#add_linking_flags_dynamic(Test)
#add_linking_flags_static(Test)
