add_library(Statistics INTERFACE
    ContingencyTable.hpp)
set_property(TARGET Statistics PROPERTY CXX_STANDARD 17)
target_include_directories(Statistics INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(Statistics-Test STATIC
    ContingencyTable.test.cpp)
target_link_libraries(Statistics-Test
    Catch2
    Linden-Load
    Linden::Genetics) #tmp

RegisterTest(ContingencyTable.test.cpp)
