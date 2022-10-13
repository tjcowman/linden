add_library(Statistics INTERFACE
    ContingencyTable.hpp)
target_include_directories(Statistics INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(Statistics-Test STATIC
    ContingencyTable.test.cpp)
target_link_libraries(Statistics-Test
    Catch2)



RegisterTest(ContingencyTable.test.cpp)

