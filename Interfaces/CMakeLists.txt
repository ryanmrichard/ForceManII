cmake_minimum_required(VERSION 3.0)
project(fmanii_api CXX)
find_package(fmanii REQUIRED)

find_package(pulsar QUIET)
if(${pulsar_FOUND})
   message(STATUS "Pulsar Framework found enabling interface.")
   add_library(fmanii_pulsar MODULE FManIIPulsarAPI.cpp
                                    creator.cpp
   )
   target_include_directories(fmanii_pulsar PRIVATE pulsar fmanii)
   target_link_libraries(fmanii_pulsar PRIVATE pulsar fmanii)
   install(TARGETS fmanii_pulsar
           LIBRARY DESTINATION ${PULSAR_MODULE_DIR}/fmanii
   )
   install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
           DESTINATION ${PULSAR_MODULE_DIR}/fmanii
           FILES_MATCHING PATTERN "*.py"
   )
endif()
