add_library(Statistics INTERFACE
    ContingencyTable.hpp)
target_include_directories(Statistics INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})