cmake_minimum_required(VERSION 3.15)

# set the project name
project(linden VERSION 1.0)



# add the executable
add_executable(linden src/main.cpp 
    src/LDForest.cpp 
    src/LDTree.cpp 
    src/Snp.cpp 
    src/TopSnpList.cpp 
    src/InputParser.cpp "include/Bitwise.h" "include/Types.h" "include/CTable.h")  

#target_compile_options(linden PRIVATE -We)

set_property(TARGET linden PROPERTY CXX_STANDARD 20)

#set_property(TARGET linden PROPERTY Wa)

#if(NOT MSCV)
   target_compile_options( linden PUBLIC -Wall -Wno-c++98-compat-pedantic  -Wno-newline-eof)
#endif()

#SET_TARGET_PROPERTIES(linden PROPERTIES COMPILE_FLAGS " -march=native" )


#set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(OpenMP REQUIRED)
target_link_libraries(linden PRIVATE OpenMP::OpenMP_CXX)

include_directories(include)


install(TARGETS linden DESTINATION bin)
