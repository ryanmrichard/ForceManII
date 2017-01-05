include_directories(${PROJECT_SOURCE_DIR})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_subdirectory(ForceFields)
add_subdirectory(InternalCoords)
add_subdirectory(ModelPotentials)
set(FMANII_SRC $<TARGET_OBJECTS:force_fields>
               $<TARGET_OBJECTS:int_coords>
               $<TARGET_OBJECTS:mod_pots>
               FManII.cpp
               ForceField.cpp
               ParameterSet.cpp
               ParseFile.cpp

)

add_library(fmanii ${FMANII_SRC})



install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" 
        DESTINATION ForceManII/include 
        FILES_MATCHING PATTERN "*.hpp"
)
