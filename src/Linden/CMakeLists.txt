
add_subdirectory(Application)
add_subdirectory(Core)
add_subdirectory(Load)

add_executable(linden
    main.cpp
)
set_property(TARGET linden PROPERTY CXX_STANDARD 17)
target_link_libraries(linden
    Linden::Application
    Linden::Core
    Linden::Load
    OpenMP::OpenMP_CXX
)

#add_executable(linden-legacy
#    main.cpp
#)
#set_property(TARGET linden PROPERTY CXX_STANDARD 17)
#target_include_directories(linden PRIVATE LindenCore LindenLoad)
#target_link_libraries(linden-legacy
#    Linden::Core::Legacy
#    Linden::Load
#    OpenMP::OpenMP_CXX
#)
#
install(TARGETS linden DESTINATION .)