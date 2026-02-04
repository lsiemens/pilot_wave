add_executable(pilot_wave
    pilot_wave.cpp)

target_sources(pilot_wave PRIVATE
    $<TARGET_OBJECTS:pilot_wave_core_obj>)

target_link_libraries(pilot_wave PRIVATE
    OpenGL::GL
    GLEW::GLEW
    glfw
    
    pilot_wave_core_obj
    pilot_wave_warnings
    pilot_wave_assets)

add_dependencies(pilot_wave pilot_wave_stage_assets)
