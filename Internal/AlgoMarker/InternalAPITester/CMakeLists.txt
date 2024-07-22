cmake_minimum_required(VERSION 3.5.0)
file(GLOB SRC_FILES
     "*.h"
     "*.cpp"
)

add_executable(InternalAPITester ${SRC_FILES})
target_link_libraries(InternalAPITester ${ADDITIONAL_LINK_FLAGS} -Wl,--start-group InfraMed Logger Mars MedAlgo MedEmbed MedIO MedMat MedPlotly MedProcessTools MedSparseMat MedSplit MedStat MedTime MedUtils QRF SerializableObject TQRF micNet -Wl,--end-group -L/nas1/Work/Libs/Boost/boost_1_85_0/installation/lib libboost_regex.so libboost_program_options.so libboost_filesystem.so libboost_system.so libboost_atomic.so /nas1/Work/SharedLibs/linux/ubuntu/${CMAKE_BUILD_TYPE}/libxgboost.so /nas1/Work/SharedLibs/linux/ubuntu/${CMAKE_BUILD_TYPE}/lib_lightgbm.so)
