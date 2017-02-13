cmake_minimum_required(VERSION 3.2)
project(fmanii-test CXX)
find_package(fmanii REQUIRED)
find_package(pulsar QUIET)
include(CTest)
enable_testing()
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


#These are the correct answers for real (large) systems
add_library(molecules SHARED
    testdata/crambin.cpp
    testdata/crambin_deriv.cpp
    testdata/peptide.cpp
    testdata/peptide_deriv.cpp
    testdata/ubiquitin.cpp
    testdata/ubiquitin_params.cpp
    testdata/ubiquitin_deriv.cpp
)
target_include_directories(molecules PRIVATE ${FMANII_ROOT})
install(TARGETS molecules DESTINATION .)

function(NEW_TEST test_name)
   add_executable(${test_name} ${test_name}.cpp) 
   target_link_libraries(${test_name} fmanii molecules)
   target_include_directories(${test_name} PRIVATE ${FMANII_ROOT})
   add_dependencies(${test_name} molecules)
   add_test(NAME ${test_name} COMMAND ${test_name})
   install(TARGETS ${test_name} DESTINATION .)
endfunction()

NEW_TEST(TestAMBER99)
NEW_TEST(TestAngle)
NEW_TEST(TestAssignParams)
NEW_TEST(TestCHARMM22)
NEW_TEST(TestDistance)
NEW_TEST(TestCoulomb)
NEW_TEST(TestFourierSeries)
NEW_TEST(TestHO)
NEW_TEST(TestLJ)
NEW_TEST(TestOPLSAA)
NEW_TEST(TestParse)
NEW_TEST(TestTorsion)
if(${pulsar_FOUND})
    include(CTestMacros)
    add_subdirectory(Interfaces)
    install(FILES ${CMAKE_BINARY_DIR}/Interfaces/CTestTestfile.cmake
            DESTINATION Interfaces)
endif()


install(FILES ${CMAKE_BINARY_DIR}/CTestTestfile.cmake DESTINATION .)

