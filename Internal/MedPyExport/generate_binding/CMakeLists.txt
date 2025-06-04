cmake_minimum_required(VERSION 3.10)
project(medpython)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/Linux/${CMAKE_BUILD_TYPE}")

set(BOOST_LIBS
    regex
    filesystem
    system
    program_options
) # for certain version, you might need to add "atomic" also

set (MEDIAL_INTERNAL_LIBS InfraMed Logger MedAlgo MedEmbed MedIO MedMat MedPlotly MedProcessTools MedSparseMat MedSplit MedStat MedTime MedUtils QRF SerializableObject TQRF micNet CommonLib
)

#set(BOOST_ROOT "$ENV{HOME}/boost-pic-install") #If you have compiled boost by yourself with -fPIC flag, please provide the path to folder with "/lib" with shared and static lib and "/include" for headers.

execute_process(COMMAND "date" "+%d-%m-%Y/%H:%M:%S"
    OUTPUT_VARIABLE DATE_TIME)
  string(REGEX REPLACE "([^/]*)/.*" "\\1"
    PROJECT_BUILD_DATE "${DATE_TIME}")
  string(REGEX REPLACE "[^/]*/([0-9:]*).*" "\\1"
    PROJECT_BUILD_TIME "${DATE_TIME}")

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} --std=c++17 -Wno-write-strings -Wuninitialized -Wno-narrowing -fopenmp -msse2 -msse3 -msse4 -march=native -Ofast --param inline-unit-growth=1000000 -fPIC -DGIT_HEAD_VERSION='\"'$(GIT_HEAD_VERSION)'\"' ")

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

# Find Boost
find_package(OpenMP REQUIRED)

find_package(Boost REQUIRED COMPONENTS ${BOOST_LIBS})  # Example components


if (BOOST_ROOT)
    set(Boost_INCLUDE_DIRS "${BOOST_ROOT}/include")
endif()
include_directories(${Boost_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../External)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../Internal)
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../External/xgboost/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../External/xgboost/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../External/xgboost/rabit/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../External/xgboost/dmlc-core/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../External/LightGBM_2.2.3/LightGBM-2.2.3/include))

# Add Internal libraries
foreach(mes_lib IN LISTS MEDIAL_INTERNAL_LIBS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../${mes_lib}/${mes_lib} ${mes_lib})
endforeach()

add_subdirectory(MedPython)
