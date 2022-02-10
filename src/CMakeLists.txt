add_subdirectory(LindenCore)

# main executable
add_executable(linden 
    linden.cpp 
) 
set_property(TARGET linden PROPERTY CXX_STANDARD 17) 
target_include_directories(linden PRIVATE LindenCore)
target_link_libraries(linden
    Linden-Core
    OpenMP::OpenMP_CXX
)

# add the original LDTree version
add_executable(lindenV0
    linden.cpp 
)  
set_property(TARGET lindenV0 PROPERTY CXX_STANDARD 17)
target_include_directories(lindenV0 PRIVATE LindenCore)
target_link_libraries(lindenV0
    Linden-Core-Old
    OpenMP::OpenMP_CXX
)

# Serializing program
add_executable(linden-serialize 
    linden-serialize.cpp
)
set_property(TARGET linden-serialize PROPERTY CXX_STANDARD 17)
target_include_directories(linden-serialize PRIVATE LindenCore)
target_link_libraries(linden-serialize
    Linden-Core
    OpenMP::OpenMP_CXX
)