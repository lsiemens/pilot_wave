find_package(Doxygen)

if(DOXYGEN_FOUND)
    set(DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/doxygen")
    set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${PROJECT_SOURCE_DIR}/README.md")
    set(DOXYGEN_INPUT_DIRS
        "${PROJECT_SOURCE_DIR}/README.md"
        "${PROJECT_SOURCE_DIR}/src")
    string(JOIN " " DOXYGEN_INPUT "${DOXYGEN_INPUT_DIRS}")

    configure_file(
        "${PROJECT_SOURCE_DIR}/doc/Doxyfile"
        "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile" @ONLY)

    doxygen_add_docs(pilot_wave_docs ALL
        "${PROJECT_SOURCE_DIR}/README.md"
        "${PROJECT_SOURCE_DIR}/src"
        CONFIG_FILE "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile"
        COMMENT "Generating Doxygen documentation")
else()
    message(WARNING "Doxygen not found")
endif()

