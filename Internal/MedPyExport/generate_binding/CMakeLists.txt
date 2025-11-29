cmake_minimum_required(VERSION 3.10)
project(medpython)
cmake_policy(SET CMP0074 NEW)
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.27")
    cmake_policy(SET CMP0144 NEW)
endif()

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

find_package(Boost REQUIRED)
if(Boost_VERSION VERSION_GREATER_EQUAL "1.69") # 1.69 is represented as 106900 in CMake
    # Remove "system" for newer versions (>= 1.69)
    # The list(REMOVE_ITEM) command is used to remove a specific item from a list variable.
    list(REMOVE_ITEM BOOST_LIBS "system")
    message(STATUS "Boost version is ${Boost_VERSION} (>= 1.69), removing 'system' from BOOST_LIBS.")
else()
    message(STATUS "Boost version is ${Boost_VERSION} (< 1.69), keeping 'system' in BOOST_LIBS.")
endif()

set (MEDIAL_INTERNAL_LIBS InfraMed Logger MedAlgo MedEmbed MedIO MedMat MedPlotly MedProcessTools MedSparseMat MedSplit MedStat MedTime MedUtils QRF SerializableObject TQRF micNet CommonLib
)

execute_process(COMMAND "date" "+%d-%m-%Y/%H:%M:%S"
    OUTPUT_VARIABLE DATE_TIME)
  string(REGEX REPLACE "([^/]*)/.*" "\\1"
    PROJECT_BUILD_DATE "${DATE_TIME}")
  string(REGEX REPLACE "[^/]*/([0-9:]*).*" "\\1"
    PROJECT_BUILD_TIME "${DATE_TIME}")

set(BUILD_STATIC_LIB ON CACHE BOOL "Build XGBoost as a static library" FORCE)
if (NOT WIN32)
    set(COMPILER_FLAGS "--std=c++17 -fopenmp -msse2 -msse3 -msse4 -march=x86-64 --param inline-unit-growth=1000000 -fPIC -Wno-write-strings -Wuninitialized -Wno-narrowing -DGIT_HEAD_VERSION='\"'$(GIT_HEAD_VERSION)'\"' ")
else()
	set(COMPILER_FLAGS "/std:c++17 /openmp /EHsc /MD /permissive- /Zc:strictStrings- /Zc:__cplusplus")
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

# Try to find static Boost libraries first
if(DEFINED ENV{BOOST_DISABLE_STATIC})
    message("BOOST_DISABLE_STATIC was found - disabling static Boost libraries")
    set(USE_BOOST_SHARED_LIBS ON)
    set(Boost_USE_STATIC_LIBS OFF)

    find_package(Boost REQUIRED COMPONENTS ${BOOST_LIBS})
else()
    set(Boost_USE_STATIC_LIBS ON)
    message("Trying to use Boost static libraries. If not compiled with -fPIC it will fail. If it is only for your local machine. You might want to set BOOST_DISABLE_STATIC environment variable")

    find_package(Boost COMPONENTS ${BOOST_LIBS})

    if(NOT Boost_FOUND)
        message(WARNING "Could not find static Boost libraries. Falling back to shared libraries.")
        set(USE_BOOST_SHARED_LIBS ON)
        set(Boost_USE_STATIC_LIBS OFF)
        # Unset Boost variables to allow a clean re-search
        unset(Boost_LIBRARIES CACHE)
        unset(Boost_FOUND CACHE)
        find_package(Boost REQUIRED COMPONENTS ${BOOST_LIBS})
    else()
        set(USE_BOOST_SHARED_LIBS OFF)
        message(INFO " Using static Boost libraries.")
    endif()
endif()




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
