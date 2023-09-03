add_library(Linden-Load STATIC
    InputParser.cpp
    InputParser.hpp
)
target_link_libraries(Linden-Load PUBLIC
    Genetics)
target_include_directories(Linden-Load PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET Linden-Load PROPERTY CXX_STANDARD 17)
add_library(Linden::Load ALIAS Linden-Load)