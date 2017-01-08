add_library(force_fields OBJECT
    amber99.cpp
    charmm22.cpp
    oplsaa.cpp
)
set_source_files_properties(amber99.cpp PROPERTIES COMPILE_FLAGS -O0)
set_source_files_properties(charmm22.cpp PROPERTIES COMPILE_FLAGS -O0)
set_source_files_properties(oplsaa.cpp PROPERTIES COMPILE_FLAGS -O0)
