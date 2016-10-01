cmake_minimum_required(VERSION 3.0)
project(ForceManII CXX)
option(BUILD_SHARED_LIBS "Should ForceManII library be shared?" ON)

#Requires C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_subdirectory(ForceManII)
add_subdirectory(tests)

