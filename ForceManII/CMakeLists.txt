set(FMANII_SRC Distance.cpp
               FManII.cpp
               HarmonicOscillator.cpp
)

add_library(fmanii ${FMANII_SRC})
target_include_directories(fmanii PUBLIC ${PROJECT_SOURCE_DIR})


install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" 
        DESTINATION ForceManII/include 
        FILES_MATCHING PATTERN "*.hpp"
)
