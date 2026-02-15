add_library(pilot_wave_quantum_obj OBJECT
    qstate.cpp square_well.cpp)

target_include_directories(pilot_wave_quantum_obj PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(pilot_wave_quantum_obj PRIVATE pilot_wave_warnings)
