if(TARGET simde)
    return()
endif()

message(STATUS "Third-party: creating target 'simde'")

include(FetchContent)
FetchContent_Declare(
    simde
    GIT_REPOSITORY https://github.com/simd-everywhere/simde.git
    GIT_TAG v0.7.2
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(simde)


add_library(simde INTERFACE)
target_include_directories(simde INTERFACE "${simde_SOURCE_DIR}")