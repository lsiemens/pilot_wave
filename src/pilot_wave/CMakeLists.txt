add_library(app STATIC)

target_sources(app PRIVATE
    qparticles.cpp controls.cpp)

target_include_directories(app PUBLIC "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(app PUBLIC core quantum
                          PRIVATE pilot_wave_warnings)


add_executable(pilot_wave)

target_sources(pilot_wave PRIVATE pilot_wave.cpp)

target_link_libraries(pilot_wave PRIVATE
    OpenGL::OpenGL
    GLEW::GLEW
    glfw
    glm::glm
    
    app
    pilot_wave_warnings)

set(ASSET_FILES
    "${CMAKE_SOURCE_DIR}/assets/vertex.glsl"
    "${CMAKE_SOURCE_DIR}/assets/fragment.glsl")
set(ASSET_DIR "${CMAKE_CURRENT_BINARY_DIR}/assets")

add_custom_command(TARGET pilot_wave POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory "${ASSET_DIR}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${ASSET_FILES} # This argument should not be in quotes so the list expands to multiple arguments
        "${ASSET_DIR}")
