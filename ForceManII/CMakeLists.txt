add_subdirectory(ForceFields)
add_subdirectory(InternalCoords)
add_subdirectory(ModelPotentials)
set(FMANII_SRC ${FMANII_SRC}
               FManII.cpp
               ParseFile.cpp
)

add_library(fmanii ${FMANII_SRC})
target_include_directories(fmanii PUBLIC ${PROJECT_SOURCE_DIR})


install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" 
        DESTINATION ForceManII/include 
        FILES_MATCHING PATTERN "*.hpp"
)
