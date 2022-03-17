cmake_minimum_required(VERSION 3.10)

project(indirectPredicates)

option(BUILD_EXECUTABLE "Build Indirect Predicates Test Executable" ON)
option(BUILD_LIBRARY    "Build Indirect Predicates Static Library"  ON)

if(BUILD_EXECUTABLE)
    message("Build Indirect Predicates Test Executable")
    add_executable(${PROJECT_NAME}-test
                   implicit_point.h
                   implicit_point.cpp
                   numerics.h
                   numerics.cpp
                   predicates/indirect_predicates.h
                   predicates/indirect_predicates.cpp
                   predicates/hand_optimized_predicates.cpp)
endif()

if(BUILD_LIBRARY)
    message("Build Indirect Predicates Static Library")
    add_library(${PROJECT_NAME}
                implicit_point.h
                implicit_point.cpp
                numerics.h
                numerics.cpp
                predicates/indirect_predicates.h
                predicates/indirect_predicates.cpp
                predicates/hand_optimized_predicates.cpp)
endif()

target_include_directories(${PROJECT_NAME} PUBLIC .)
target_include_directories(${PROJECT_NAME} PUBLIC predicates)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_11)

# Compiler-specific options to grant IEEE 754 compliance
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(${PROJECT_NAME} PUBLIC -frounding-math)
    find_package(Threads)
    target_link_libraries(${PROJECT_NAME} ${CMAKE_THREAD_LIBS_INIT})
    #target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    target_compile_options(${PROJECT_NAME} PUBLIC "/fp:strict")
    target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")
    target_compile_options(${PROJECT_NAME} PUBLIC "/D _CRT_SECURE_NO_WARNINGS")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR
       "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    #target_compile_options(${PROJECT_NAME} PUBLIC -O0)
endif()
