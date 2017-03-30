cmake_minimum_required(VERSION 3.0)
project(ForceManII-SuperBuild CXX)

include(ExternalProject)

option(BUILD_SHARED_LIBS "Should ForceManII library be shared?" ON)

#Requires C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(STAGE_DIR            "${CMAKE_BINARY_DIR}/stage")
set(STAGE_INSTALL_PREFIX "${STAGE_DIR}${CMAKE_INSTALL_PREFIX}")
set(FMANII_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")
ExternalProject_Add(fmanii
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ForceManII
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DFMANII_ROOT=${FMANII_ROOT}
    BUILD_ALWAYS 1
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install DESTDIR=${STAGE_DIR}
    CMAKE_CACHE_ARGS -DCMAKE_INSTALL_RPATH:LIST=${CMAKE_INSTALL_RPATH}
                     -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                     -DCMAKE_PREFIX_PATH:LIST=${CMAKE_PREFIX_PATH}
)
list(APPEND FMANII_PREFIX_PATH ${STAGE_INSTALL_PREFIX} ${CMAKE_PREFIX_PATH})

find_package(pulsar QUIET)
if(${pulsar_FOUND})
   message(STATUS "Pulsar Framework found enabling interface.")
    ExternalProject_Add(pulsar_api
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Interfaces
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                   -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        BUILD_ALWAYS 1
        INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install DESTDIR=${STAGE_DIR}
        CMAKE_CACHE_ARGS -DCMAKE_PREFIX_PATH:LIST=${FMANII_PREFIX_PATH}
                         -DCMAKE_INSTALL_RPATH:LIST=${CMAKE_INSTALL_RPATH}
                         -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    )
    add_dependencies(pulsar_api fmanii)
else()
    add_library(pulsar_api INTERFACE)
endif()


ExternalProject_Add(fmanii_test
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/test_stage
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DFMANII_ROOT=${FMANII_ROOT}
               -DSTAGE_DIR=${STAGE_DIR}
    BUILD_ALWAYS 1
    CMAKE_CACHE_ARGS -DCMAKE_PREFIX_PATH:LIST=${FMANII_PREFIX_PATH}
                     -DCMAKE_INSTALL_RPATH:LIST=${CMAKE_INSTALL_RPATH}
                     -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
)
add_dependencies(fmanii_test fmanii pulsar_api)

file(WRITE ${CMAKE_BINARY_DIR}/CTestTestfile.cmake "subdirs(test_stage)")

install(DIRECTORY ${CMAKE_BINARY_DIR}/stage/${CMAKE_INSTALL_PREFIX}/
        DESTINATION ${CMAKE_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS)


