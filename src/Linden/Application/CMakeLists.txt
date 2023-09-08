
add_library(Linden-Application STATIC
    Ping.hpp
    Session.cpp
    Session.hpp
    Listener.cpp
    Listener.hpp
    Platform.cpp
    Platform.hpp)
target_include_directories(Linden-Application PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET Linden-Application PROPERTY CXX_STANDARD 17)
add_library(Linden::Application ALIAS Linden-Application)