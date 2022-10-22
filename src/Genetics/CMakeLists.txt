add_library(Genetics INTERFACE
    GenotypeMatrix.hpp
    Types.hpp
    Location.hpp
    Locus.hpp)
target_include_directories(Genetics INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(Linden::Genetics ALIAS Genetics)