### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# project
cmake_minimum_required(VERSION 3.15...4.0)
project(RidgeSurfaceFinder VERSION 1.0
        DESCRIPTION "Ridge Surface Extraction via Seed Points"
        LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
configure_file(Config.h.in Config.h)

# definitions
add_definitions(-D__DATAPATH_RAW__="${CMAKE_CURRENT_BINARY_DIR}/data")

# create library
add_subdirectory(includes/boost_heap)
add_subdirectory(includes/boost_hash)
add_subdirectory(includes/fmt)

# add_subdirectory(src)
set(SOURCES src/utils/Vec.cpp src/rsf/CubicalRidgeSurfaceFinder.cpp src/io/npy.cpp src/rsf/ChangeGraph.cpp)
add_library(rsflib STATIC ${SOURCES})
target_include_directories(rsflib PUBLIC src)
target_link_libraries(rsflib PUBLIC fmt::fmt-header-only)
target_include_directories(rsflib PUBLIC includes/parallel-hashmap)
target_include_directories(rsflib PRIVATE includes/eigen)

# create main exec
add_executable(rsf main.cxx)
target_link_libraries(rsf PUBLIC rsflib)

# create examples
add_executable(examples examples.cxx)
target_link_libraries(examples PUBLIC rsflib)

# create tests
add_executable(tests test.cxx)
add_subdirectory(includes/catch2)
target_link_libraries(tests PUBLIC rsflib)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)

# create C API
add_library(rsfapi SHARED api.cpp)
target_link_libraries(rsfapi PUBLIC rsflib)

# copy test data
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})