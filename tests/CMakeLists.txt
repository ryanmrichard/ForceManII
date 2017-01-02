set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include(CTest)
enable_testing()

#These are the correct answers
add_library(molecules SHARED testdata/ubiquitin.cpp
                             testdata/ubiquitin_intcoords.cpp
                             testdata/ubiquitin_params.cpp
)
target_include_directories(molecules PRIVATE ${PROJECT_SOURCE_DIR})
add_dependencies(molecules fmanii)

function(NEW_TEST test_name)
   add_executable(${test_name} ${test_name}.cpp) 
   target_link_libraries(${test_name} fmanii molecules)
   target_include_directories(${test_name} PRIVATE ${PROJECT_SOURCE_DIR})
   add_dependencies(${test_name} fmanii molecules)
   add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

NEW_TEST(TestAssignParams)
NEW_TEST(TestCoulomb)
NEW_TEST(TestFourierSeries)
NEW_TEST(TestHO)
NEW_TEST(TestIntCoords)
NEW_TEST(TestLJ)
NEW_TEST(TestParse)

