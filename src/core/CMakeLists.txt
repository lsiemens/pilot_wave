add_library(core STATIC)

target_sources(core PRIVATE
    looplog.cpp frame_timer.cpp
    model.cpp camera.cpp rng.cpp
    object.cpp particles.cpp shaders.cpp
    geometry.cpp math_util.cpp)

target_include_directories(core PUBLIC "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(core PRIVATE pilot_wave_warnings)
