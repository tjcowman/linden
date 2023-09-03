cmake_minimum_required(VERSION 3.24)

# set the project name
project(linden VERSION 1.0)

#if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
#endif()

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)

# Function to add the cpp file to the global list of tests
include(cmake/RegisterTest.cmake)

add_subdirectory(deps)
add_subdirectory(src)
add_subdirectory(example)
add_subdirectory(test)

#install(TARGETS linden lindenV0 DESTINATION .)
