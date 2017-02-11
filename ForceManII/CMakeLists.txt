cmake_minimum_required(VERSION 3.0)
project(fmanii CXX)
include_directories(${FMANII_ROOT})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_subdirectory(ForceFields)
add_subdirectory(InternalCoords)
add_subdirectory(ModelPotentials)
set(FMANII_SRC $<TARGET_OBJECTS:force_fields>
               $<TARGET_OBJECTS:int_coords>
               $<TARGET_OBJECTS:mod_pots>
               FManII.cpp
               FFTerm.cpp
               ForceField.cpp
               ParameterSet.cpp
               ParseFile.cpp
)
add_library(fmanii ${FMANII_SRC})

install(TARGETS fmanii
        DESTINATION lib)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" 
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp"
)
configure_file("fmaniiConfig.cmake.in" fmaniiConfig.cmake @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/fmaniiConfig.cmake
        DESTINATION share/cmake/fmanii)
