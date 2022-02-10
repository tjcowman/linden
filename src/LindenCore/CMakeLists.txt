

add_library(Linden-Core STATIC
    CommonStructs.h
    CTable.cpp
    CTable.h
    Bitwise.cpp
    Bitwise.h
    Graph.cpp
    Graph.h
    InputParser.h
    LDForest.cpp
    LDForest.h
    LDTree2.cpp
    LDTree2.h
    Snp.cpp
    Snp.h
    SnpSet.cpp
    SnpSet.h
    Types.h
    TopSnpList.cpp
    TopSnpList.h
)
set_property(TARGET Linden-Core PROPERTY CXX_STANDARD 17)
target_link_libraries(Linden-Core 
    OpenMP::OpenMP_CXX
 ) 


add_library(Linden-Core-Old STATIC
    CommonStructs.h
    CTable.cpp
    CTable.h
    Bitwise.cpp
    Bitwise.h
    Graph.cpp
    Graph.h
    InputParser.h
    LDForest.cpp
    LDForest.h
    LDTree.cpp
    LDTree.h
    Snp.cpp
    Snp.h
    SnpSet.cpp
    SnpSet.h
    Types.h
    TopSnpList.cpp
    TopSnpList.h
)
target_compile_definitions(Linden-Core-Old PUBLIC -DoldTree)
set_property(TARGET Linden-Core-Old PROPERTY CXX_STANDARD 17) 
target_link_libraries(Linden-Core-Old  
    OpenMP::OpenMP_CXX
) 