set(FMANII_SRC Angle.cpp
               Distance.cpp
               FManII.cpp
               FourierSeries.cpp
               HarmonicOscillator.cpp
               Torsion.cpp
)

add_library(fmanii ${FMANII_SRC})
target_include_directories(fmanii PUBLIC ${PROJECT_SOURCE_DIR})


install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" 
        DESTINATION ForceManII/include 
        FILES_MATCHING PATTERN "*.hpp"
)
