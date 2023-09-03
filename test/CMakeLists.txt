
get_property(reg_tests GLOBAL PROPERTY reg_tests)

add_executable(linden-tests
    main.cpp
    ${reg_tests})

set_property(TARGET linden-tests PROPERTY CXX_STANDARD 17)

target_link_libraries(linden-tests PUBLIC
    Catch2
    Linden-Load)
    
install(TARGETS linden-tests DESTINATION test)
