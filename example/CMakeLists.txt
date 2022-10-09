
configure_file(cases.in cases)
configure_file(controls.in controls)
configure_file(loci.in loci)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/cases
    ${CMAKE_CURRENT_BINARY_DIR}/controls
    ${CMAKE_CURRENT_BINARY_DIR}/loci
    DESTINATION /)