cmake_minimum_required(VERSION 3.2)
project(fmanii-test CXX)
find_package(fmanii REQUIRED)
find_package(pulsar QUIET)
include(CTest)
enable_testing()

#These are the correct answers
add_library(molecules SHARED
    testdata/crambin.cpp
    testdata/crambin_intcoords.cpp
    testdata/peptide.cpp
    testdata/peptide_intcoords.cpp
    testdata/ubiquitin.cpp
    testdata/ubiquitin_intcoords.cpp
    testdata/ubiquitin_params.cpp
)
target_include_directories(molecules PRIVATE ${FMANII_ROOT})

function(NEW_TEST test_name)
   add_executable(${test_name} ${test_name}.cpp) 
   target_link_libraries(${test_name} fmanii molecules)
   target_include_directories(${test_name} PRIVATE ${FMANII_ROOT})
   add_dependencies(${test_name} molecules)
   add_test(NAME ${test_name} COMMAND ${test_name})
   install(TARGETS ${test_name} DESTINATION .)
endfunction()

NEW_TEST(TestAssignParams)
NEW_TEST(TestDeriv)
NEW_TEST(TestCoulomb)
NEW_TEST(TestFourierSeries)
NEW_TEST(TestHO)
NEW_TEST(TestIntCoords)
NEW_TEST(TestLJ)
NEW_TEST(TestParse)
if(${pulsar_FOUND})
    include(CTestMacros)
    pulsar_py_test(. TestPulsarAPI ${STAGE_DIR}${PULSAR_MODULE_DIR})
endif()


install(FILES ${CMAKE_BINARY_DIR}/CTestTestfile.cmake DESTINATION .)

