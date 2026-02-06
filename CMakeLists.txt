cmake_minimum_required(VERSION 3.19)
project(pilot_wave
    VERSION 1.0
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# setup OpenGL GLEW and GLFW
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(glfw3 CONFIG)

option(ENABLE_WARNINGS "Enable compiler warnings" ON)
option(STRICT_WARNINGS "Enable strict compiler warnings" ON)
option(ERROR_ON_WARNING "Treat warnings to errors" OFF)
option(BUILD_DOC "Build documentation" ON)

add_library(pilot_wave_warnings INTERFACE)
if(ENABLE_WARNINGS)
    target_compile_options(pilot_wave_warnings INTERFACE
        $<$<CXX_COMPILER_ID:GNU,Clang>: -Wall -Wextra -Wpedantic -Wshadow -Wformat=2 -Wundef>
        $<$<AND:$<BOOL:${STRICT_WARNINGS}>,$<CXX_COMPILER_ID:GNU,Clang>>: -Wconversion>
        $<$<CXX_COMPILER_ID:Clang>: -Wdocumentation>

        $<$<AND:$<BOOL:${ERROR_ON_WARNING}>,$<CXX_COMPILER_ID:GNU,Clang>>: -Werror>)
endif()

if(BUILD_DOC)
    add_subdirectory(doc)
endif()

add_subdirectory(src)
add_subdirectory(assets)
