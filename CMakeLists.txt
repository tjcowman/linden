cmake_minimum_required(VERSION 3.20)

# set the project name
project(linden VERSION 1.0)

find_package(OpenMP REQUIRED)
find_package(Catch2 3 REQUIRED)

add_subdirectory(src)



#add_executable(linden-tests test/main.cpp
#    include/LDTree2.h
#    src/LDTree2.cpp 
#)
#set_property(TARGET linden-tests PROPERTY CXX_STANDARD 17)
#target_link_libraries(linden-tests
#    Linden-Base
#   OpenMP::OpenMP_CXX
#    Catch2::Catch2WithMain
#)

install(TARGETS linden DESTINATION bin)
