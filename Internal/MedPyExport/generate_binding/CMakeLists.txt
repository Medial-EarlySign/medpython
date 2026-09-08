cmake_minimum_required(VERSION 3.10)
project(medpython)
cmake_policy(SET CMP0074 NEW)
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.27")
    cmake_policy(SET CMP0144 NEW)
endif()
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.13")
    cmake_policy(SET CMP0077 NEW)
endif()

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
            SO_COMPILATION
        )
    endif()
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find Boost
find_package(OpenMP REQUIRED)

find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
include(FetchContent)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/Linux/${CMAKE_BUILD_TYPE}")
set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR})

if(DEFINED ENV{BOOST_ROOT})
    message(STATUS "BOOST_ROOT is set to: $ENV{BOOST_ROOT}")
    set(BOOST_ROOT "$ENV{BOOST_ROOT}")
    set(Boost_INCLUDE_DIRS "${BOOST_ROOT}/include")
    include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
else()
    find_package(Boost QUIET)
    if (NOT Boost_FOUND)
        message(STATUS "BOOST_ROOT environment variable is not set. Please set it to the Boost installation path, if you want to use a specific boost.")
        message(STATUS "Boost not found. Downloading headers...")
        
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
        if(NOT Boost_INCLUDE_DIRS AND TARGET Boost::boost)
            get_target_property(Boost_INCLUDE_DIRS Boost::boost INTERFACE_INCLUDE_DIRECTORIES)
        endif()
        message(STATUS "Boost found. Using system headers at: ${Boost_INCLUDE_DIRS}")
        include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
    endif()
endif()

if(POLICY CMP0091)
    cmake_policy(SET CMP0091 NEW)
endif()
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL" CACHE STRING "Runtime library" FORCE)

## Eigen
find_package(Eigen3 QUIET)

if(Eigen3_FOUND)
    if(NOT EIGEN3_INCLUDE_DIR AND TARGET Eigen3::Eigen)
        get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    endif()
    message(STATUS "Eigen3 found. Using system headers at: ${EIGEN3_INCLUDE_DIR}")
    include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
else()
    message(STATUS "Eigen3 not found. Downloading headers from source...")
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
    if(NOT nlohmann_json_INCLUDE_DIR AND TARGET nlohmann_json::nlohmann_json)
        get_target_property(nlohmann_json_INCLUDE_DIR nlohmann_json::nlohmann_json INTERFACE_INCLUDE_DIRECTORIES)
    endif()
    message(STATUS "nlohmann_json found. Using system headers at: ${nlohmann_json_INCLUDE_DIR}")
    include_directories(SYSTEM ${nlohmann_json_INCLUDE_DIR})
else ()
    # https://github.com/nlohmann/json/releases/download/v3.11.2/json.hpp include.zip
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

set(XGBOOST_STATIC_FALLBACK OFF)
find_package(xgboost QUIET)

if(NOT xgboost_FOUND)
    find_path(XGBOOST_INCLUDE_DIR NAMES xgboost/c_api.h)
    find_library(XGBOOST_LIB NAMES xgboost libxgboost)

    if(XGBOOST_INCLUDE_DIR AND XGBOOST_LIB)
        message(STATUS "XGBoost CMake config missing, but raw binaries found!")
        add_library(xgboost::xgboost UNKNOWN IMPORTED)
        set_target_properties(xgboost::xgboost PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${XGBOOST_INCLUDE_DIR}"
            IMPORTED_LOCATION "${XGBOOST_LIB}"
        )
        set(xgboost_FOUND TRUE)
    endif()
endif()

if(xgboost_FOUND)
    message(STATUS "XGBoost found via find_package. Using shared library.")
    set(XGBOOST_TARGET xgboost::xgboost)
else()
    message(STATUS "XGBoost not found. Falling back to static FetchContent...")
    set(XGBOOST_STATIC_FALLBACK ON)
    FetchContent_Declare(
        xgboost
        URL      https://github.com/dmlc/xgboost/releases/download/v${XGBOOST_VERSION}/xgboost-${XGBOOST_VERSION}.tar.gz
        DOWNLOAD_EXTRACT_TIMESTAMP ON
        SOURCE_SUBDIR ${XGBOOST_TARGET_DIR}
    )

    set(BUILD_STATIC_LIB ON CACHE BOOL "Build XGBoost as a static library" FORCE)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static libs" FORCE)

    if(WIN32)
        # --- FORCE XGBOOST & DMLC TO USE /MD ---
        set(FORCE_SHARED_CRT ON CACHE BOOL "Disable XGBoost CLI executable" FORCE)
    endif()
    FetchContent_MakeAvailable(xgboost)

    include_directories(SYSTEM ${xgboost_SOURCE_DIR}/cpp_src/include)          
    include_directories(SYSTEM ${xgboost_SOURCE_DIR}/cpp_src/rabit/include)
    include_directories(SYSTEM ${xgboost_SOURCE_DIR}/cpp_src/dmlc-core/include)

    set(XGBOOST_TARGET xgboost)
endif()

# LightGBM
set(LIGHTGBM_STATIC_FALLBACK OFF)
find_package(LightGBM QUIET)

if(NOT LightGBM_FOUND)
    find_path(LIGHTGBM_INCLUDE_DIR NAMES LightGBM/c_api.h)
    find_library(LIGHTGBM_LIB NAMES lightgbm lib_lightgbm _lightgbm)

    if(LIGHTGBM_INCLUDE_DIR AND LIGHTGBM_LIB)
        message(STATUS "LightGBM CMake config missing, but raw binaries found!")
        add_library(LightGBM::lightgbm UNKNOWN IMPORTED)
        set_target_properties(LightGBM::lightgbm PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${LIGHTGBM_INCLUDE_DIR}"
            IMPORTED_LOCATION "${LIGHTGBM_LIB}"
        )
        set(LightGBM_FOUND TRUE)
    endif()
endif()

if(LightGBM_FOUND)
    message(STATUS "LightGBM found via find_package. Using shared library.")
    set(LIGHTGBM_TARGET LightGBM::lightgbm)
else()
    message(STATUS "LightGBM not found. Falling back to static FetchContent...")
    set(LIGHTGBM_STATIC_FALLBACK ON)
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
    
    
    set(LIGHTGBM_TARGET _lightgbm)
endif()

## End

set (MEDIAL_INTERNAL_LIBS InfraMed Logger MedAlgo MedEmbed MedIO MedMat MedPlotly MedProcessTools MedSparseMat MedSplit MedStat MedTime MedUtils QRF SerializableObject TQRF micNet CommonLib
)

if (NOT WIN32)
    # 1. Base flags common to all Unix/Linux/Mac systems
    set(COMMON_FLAGS "-fPIC -Wno-write-strings -Wuninitialized -Wno-narrowing -DMES_LIBRARY=1 -DGIT_HEAD_VERSION='\"'$(GIT_HEAD_VERSION)'\"'")

    if(APPLE)
        set(COMPILER_FLAGS "${COMMON_FLAGS} -Xpreprocessor -fopenmp -D_LIBCPP_DISABLE_AVAILABILITY")
    else()
        # Linux specific base flags
        if(CMAKE_SIZEOF_VOID_P EQUAL 4)
            message(STATUS "32-bit architecture detected. Restricting compiler memory usage...")
            set(COMPILER_FLAGS "${COMMON_FLAGS} -g0 -ftrack-macro-expansion=0 -fopenmp")
        else()
            set(COMPILER_FLAGS "${COMMON_FLAGS} -fopenmp --param inline-unit-growth=1000000")
        endif()


        # 2. Check Architecture
        # Matches aarch64 (Linux) or arm64 (Apple Silicon)
        if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64)$")
            message(STATUS "Target architecture: x86_64")
            set(COMPILER_FLAGS "${COMPILER_FLAGS} -msse2 -msse3 -msse4 -march=x86-64")
        elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(i686)$")
            message(STATUS "Detected x86 (32-bit) architecture - i686")
            set(COMPILER_FLAGS "${COMPILER_FLAGS} -msse2 -march=i686")
        elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(i386|x86)$")
            message(STATUS "Detected x86 (32-bit) architecture")
            set(COMPILER_FLAGS "${COMPILER_FLAGS} -march=i386")
        else()
            message(STATUS "Using Generic Flags To Target architecture: ${CMAKE_SYSTEM_PROCESSOR}")
        endif()
    endif()
else()
    # Windows
    set(COMPILER_FLAGS "/openmp /EHsc /MD /permissive- /Zc:strictStrings- /Zc:__cplusplus /utf-8 /DMES_LIBRARY=1")
    add_compile_definitions(GIT_HEAD_VERSION="$ENV{GIT_HEAD_VERSION}")
    add_compile_definitions(BOOST_ALL_NO_LIB)
endif()
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${COMPILER_FLAGS} ")
# Filter out definition of NDEBUG from the default configuration flags for Release
foreach (language CXX C)
	set(VAR_TO_MODIFY "CMAKE_${language}_FLAGS_RELEASE")
	string(REGEX REPLACE "(^| )[/-]D *NDEBUG($| )"
						 " "
						 replacement
						 "${${VAR_TO_MODIFY}}"
		  )
	message("Original (${VAR_TO_MODIFY}) is ${${VAR_TO_MODIFY}} replacement is ${replacement}")
	set(${VAR_TO_MODIFY} "${replacement}" CACHE STRING "Default flags for Release configuration" FORCE)
endforeach()


set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -march=native -Og")

set(CMAKE_BUILD_TYPE "Release")


set(PYTHON_INCLUDE_DIR ${Python3_INCLUDE_DIRS})


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../Internal)
include_directories(${PROJECT_SOURCE_DIR})

# Add Internal libraries
foreach(mes_lib IN LISTS MEDIAL_INTERNAL_LIBS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../${mes_lib}/${mes_lib} ${mes_lib})
endforeach()

include_directories(../MedPyExport)

add_subdirectory(MedPython)
