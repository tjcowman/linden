cmake_minimum_required(VERSION 3.20)

# set the project name
project(linden VERSION 1.0)

#if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
#endif()

find_package(OpenMP REQUIRED)

# Function to add the cpp file to the global list of tests
include(cmake/RegisterTest.cmake)

add_subdirectory(deps)
add_subdirectory(src)
add_subdirectory(example)



add_subdirectory(test)





set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/install)

#install(TARGETS linden DESTINATION bin)
install(TARGETS linden lindenV0 DESTINATION ${CMAKE_INTALL_PREFIX}/)
