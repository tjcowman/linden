cmake_minimum_required(VERSION 3.20)

# set the project name
project(linden VERSION 1.0)

find_package(OpenMP REQUIRED)
#find_package(Catch2 REQUIRED)

add_subdirectory(deps)
add_subdirectory(src)
add_subdirectory(example)
add_subdirectory(test)





set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/install)

#install(TARGETS linden DESTINATION bin)
install(TARGETS linden lindenV0 DESTINATION ${CMAKE_INTALL_PREFIX}/)
