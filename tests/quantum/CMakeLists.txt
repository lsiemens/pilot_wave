set(TESTS qstate_test square_well_test harmonic_oscillator_test)

foreach (test_name IN LISTS TESTS)
    add_executable(${test_name} ${test_name}.cpp)

    target_link_libraries(${test_name} PRIVATE
        glm::glm
        quantum
        testing
        GTest::gtest_main
        pilot_wave_warnings)

    gtest_discover_tests(${test_name})
endforeach()
