
function(RegisterTest SOURCE)
    get_property(tmp GLOBAL PROPERTY reg_tests)
    set_property(GLOBAL PROPERTY reg_tests 
        ${tmp} ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE})
endfunction()