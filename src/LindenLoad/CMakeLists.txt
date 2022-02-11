add_library(Linden-Load STATIC
    InputParser.cpp
    InputParser.hpp
    Serializers.hpp
)
set_property(TARGET Linden-Load PROPERTY CXX_STANDARD 17)