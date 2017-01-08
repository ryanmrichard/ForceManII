cmake_minimum_required(VERSION 3.0)
project(ForceManII CXX)
option(BUILD_SHARED_LIBS "Should ForceManII library be shared?" ON)

#Requires C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

add_subdirectory(ForceManII)
add_subdirectory(tests)

#This file will allow us to run ctest in the top-level build dir
#Basically it just defers to the actual top-level CTestTestfile.cmake in the
#build directory for this project
file(WRITE ${CMAKE_BINARY_DIR}/CTestTestfile.cmake "subdirs(tests)")


