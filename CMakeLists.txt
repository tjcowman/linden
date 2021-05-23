cmake_minimum_required(VERSION 3.10)

# set the project name
project(linden VERSION 1.0)



# add the executable
add_executable(linden src/main.cpp 
    src/LDForest.cpp 
    src/LDTree.cpp 
    src/Snp.cpp 
    src/TopSnpList.cpp 
   # src/TsvParser.cpp 
    src/InputParser.cpp)  


set_property(TARGET linden PROPERTY CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(OpenMP REQUIRED)
target_link_libraries(linden PRIVATE OpenMP::OpenMP_CXX)

include_directories(include)


install(TARGETS linden DESTINATION bin)
