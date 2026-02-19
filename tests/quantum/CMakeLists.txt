set(TESTS square_well_test harmonic_oscillator_test)

foreach (test_name IN LISTS TESTS)
    add_executable(${test_name} ${test_name}.cpp)

    target_link_libraries(${test_name} PRIVATE
        glm::glm
        quantum
        testing
        pilot_wave_warnings)

    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
