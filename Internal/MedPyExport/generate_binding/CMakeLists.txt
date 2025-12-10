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

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/Linux/${CMAKE_BUILD_TYPE}")
set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR})

if(DEFINED ENV{BOOST_ROOT})
    message(STATUS "BOOST_ROOT is set to: $ENV{BOOST_ROOT}")
    set(BOOST_ROOT "$ENV{BOOST_ROOT}")
else()
    message(STATUS "BOOST_ROOT environment variable is not set. Please set it to the Boost installation path.")
    # set(BOOST_ROOT "$ENV{HOME}/Documents/MES/Boost") #If you have compiled boost by yourself with -fPIC flag, please provide the path to folder with "/lib" with shared and static lib and "/include" for headers.
endif()

set(BOOST_LIBS
    regex
    filesystem
    system
    program_options
) # for certain version, you might need to add "atomic" also, in some other cases remove system

find_package(Boost QUIET)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "Could not locate Boost installation to determine version.")
endif()

if(Boost_VERSION_STRING VERSION_GREATER_EQUAL "1.69") # 1.69 is represented as 106900 in CMake
    # Remove "system" for newer versions (>= 1.69)
    # The list(REMOVE_ITEM) command is used to remove a specific item from a list variable.
    list(REMOVE_ITEM BOOST_LIBS "system")
    message(STATUS "Boost version is ${Boost_VERSION_STRING} (>= 1.69), removing 'system' from BOOST_LIBS.")
else()
    message(STATUS "Boost version is ${Boost_VERSION_STRING} (< 1.69), keeping 'system' in BOOST_LIBS.")
endif()
if (WIN32)
    list(APPEND BOOST_LIBS "serialization")
endif()

# Try to find static Boost libraries first
if(DEFINED ENV{BOOST_DISABLE_STATIC})
    message("BOOST_DISABLE_STATIC was found - disabling static Boost libraries")
    set(USE_BOOST_SHARED_LIBS ON)
    set(Boost_USE_STATIC_LIBS OFF)

    find_package(Boost REQUIRED COMPONENTS ${BOOST_LIBS})
    if (WIN32)
        add_definitions(-DBOOST_ALL_DYN_LINK) # Dynamic
    endif()
else()
    # Attempt 1: Try to use static Boost libraries first
    set(Boost_USE_STATIC_LIBS ON)
    message("Trying to use Boost static libraries. If not compiled with -fPIC it will fail. If it is only for your local machine. You might want to set BOOST_DISABLE_STATIC environment variable")

    find_package(Boost QUIET COMPONENTS ${BOOST_LIBS})

    if(NOT Boost_FOUND)
        message(WARNING "Could not find static Boost libraries. Falling back to shared libraries.")

        set(USE_BOOST_SHARED_LIBS ON)
        set(Boost_USE_STATIC_LIBS OFF)
        # Unset Boost variables to allow a clean re-search
        unset(Boost_LIBRARIES CACHE)
        unset(Boost_FOUND CACHE)
        if (WIN32)
            add_definitions(-DBOOST_ALL_DYN_LINK) # Dynamic
        endif()
        find_package(Boost REQUIRED COMPONENTS ${BOOST_LIBS})
    else()
        set(USE_BOOST_SHARED_LIBS OFF)
        if (WIN32)
            add_definitions(-DBOOST_ALL_NO_LIB) # Static
        endif()
        message(STATUS  "Successfully found static Boost libraries.")
    endif()
endif()

set (MEDIAL_INTERNAL_LIBS InfraMed Logger MedAlgo MedEmbed MedIO MedMat MedPlotly MedProcessTools MedSparseMat MedSplit MedStat MedTime MedUtils QRF SerializableObject TQRF micNet CommonLib
)

set(BUILD_STATIC_LIB ON CACHE BOOL "Build XGBoost as a static library" FORCE)
if (NOT WIN32)
    # 1. Base flags common to all Unix/Linux/Mac systems
    set(COMMON_FLAGS "-fPIC -Wno-write-strings -Wuninitialized -Wno-narrowing -DGIT_HEAD_VERSION='\"'$(GIT_HEAD_VERSION)'\"'")

    if(APPLE)
        set(COMPILER_FLAGS "${COMMON_FLAGS} -Xpreprocessor -fopenmp")
    else()
        # Linux specific base flags
        set(COMPILER_FLAGS "${COMMON_FLAGS} -fopenmp --param inline-unit-growth=1000000")

        # 2. Check Architecture
        # Matches aarch64 (Linux) or arm64 (Apple Silicon)
        if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64)")
            message(STATUS "Target architecture: ARM64 (aarch64)")
            # No -march=x86-64 or -msse flags here! 
            # ARM64 implies NEON support by default, so usually no extra SIMD flags needed.
        else()
            message(STATUS "Target architecture: x86_64")
            # Apply Intel/AMD specific flags
            set(COMPILER_FLAGS "${COMPILER_FLAGS} -msse2 -msse3 -msse4 -march=x86-64")
        endif()
    endif()
else()
    # Windows
    set(COMPILER_FLAGS "/openmp /EHsc /MD /permissive- /Zc:strictStrings- /Zc:__cplusplus")
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


if (BOOST_ROOT)
    set(Boost_INCLUDE_DIRS "${BOOST_ROOT}/include")
endif()
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../External)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../Internal)
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../External/xgboost/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../External/xgboost/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../External/xgboost/rabit/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../External/xgboost/dmlc-core/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../External/LightGBM_2.2.3/LightGBM-2.2.3/include)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static libs" FORCE)
if(POLICY CMP0091)
    cmake_policy(SET CMP0091 NEW)
endif()
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL" CACHE STRING "Runtime library" FORCE)

# Add Internal libraries
foreach(mes_lib IN LISTS MEDIAL_INTERNAL_LIBS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../${mes_lib}/${mes_lib} ${mes_lib})
endforeach()

include_directories(../MedPyExport)

# External Libraries
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../../External/xgboost xgboost)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../../External/LightGBM_2.2.3/LightGBM-2.2.3 _lightgbm)

add_subdirectory(MedPython)