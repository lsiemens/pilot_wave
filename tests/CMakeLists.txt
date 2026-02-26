if(FULL_TEST)
    add_compile_definitions(FULL_TEST)
endif()

include(FetchContent)

FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.17.0)

set(_OLD_BULID_TESTING ${BUILD_TESTING})
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
set(BUILD_TESTING ${_OLD_BUILD_TESTING} CACHE BOOL "" FORCE)
include(GoogleTest)

add_subdirectory(core)
add_subdirectory(quantum)
