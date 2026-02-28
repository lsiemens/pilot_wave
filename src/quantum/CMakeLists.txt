add_library(quantum STATIC)

target_sources(quantum PRIVATE
    qstate.cpp square_well.cpp harmonic_oscillator.cpp hydrogen_atom.cpp)

target_include_directories(quantum PUBLIC "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(quantum PUBLIC core
                              PRIVATE pilot_wave_warnings)
