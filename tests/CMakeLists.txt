set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include(CTest)
enable_testing()

add_executable(testHO TestHO.cpp)
target_link_libraries(testHO fmanii)
target_include_directories(testHO PRIVATE ${PROJECT_SOURCE_DIR})
add_dependencies(testHO fmanii)

add_test(NAME HOTest COMMAND testHO)