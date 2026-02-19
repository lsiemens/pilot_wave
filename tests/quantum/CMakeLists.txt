add_executable(square_well_test square_well_test.cpp)

target_link_libraries(square_well_test PRIVATE
    glm::glm
    quantum
    testing
    pilot_wave_warnings)

add_test(NAME square_well_test COMMAND square_well_test)
