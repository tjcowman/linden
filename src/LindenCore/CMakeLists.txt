
add_library(Linden-Core STATIC
    CTable.h
    Bitwise.h
    Graph.h
    LDForest.cpp
    LDForest.h
    LDTree2.cpp
    LDTree2.h
    Snp.cpp
    Snp.h
    SnpSet.cpp
    SnpSet.h
    TopSnpList.cpp
    TopSnpList.h
)
set_property(TARGET Linden-Core PROPERTY CXX_STANDARD 17)
target_link_libraries(Linden-Core 
    OpenMP::OpenMP_CXX
 )
 target_include_directories(Linden-Core PRIVATE ../LindenLoad) 


add_library(Linden-Core-Old STATIC
    CTable.h
    Bitwise.h
    Graph.h
    LDForest.cpp
    LDForest.h
    LDTree.cpp
    LDTree.h
    Snp.cpp
    Snp.h
    SnpSet.cpp
    SnpSet.h
    TopSnpList.cpp
    TopSnpList.h
)
target_compile_definitions(Linden-Core-Old PUBLIC -DoldTree)
set_property(TARGET Linden-Core-Old PROPERTY CXX_STANDARD 17) 
target_link_libraries(Linden-Core-Old  
    OpenMP::OpenMP_CXX
)
target_include_directories(Linden-Core-Old PRIVATE ../LindenLoad)  