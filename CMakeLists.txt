cmake_minimum_required(VERSION 3.10)

# set the project name
project(indirectPredicates)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# add the executable
add_library(${PROJECT_NAME}
    implicit_point.cpp
    numerics.cpp
    predicates/indirect_predicates.cpp
    predicates/hand_optimized_predicates.cpp
)

# Compiler-specific options
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME} PUBLIC -frounding-math)
# use intrinsic functions (CHECK WHAT TO DO FOR GCC !!!!!!!!)
#	target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME} PUBLIC "/fp:strict")
# use intrinsic functions
	target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")
# turn off annoying warnings
	target_compile_options(${PROJECT_NAME} PUBLIC "/D _CRT_SECURE_NO_WARNINGS")
endif()

# Public include directory
target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/predicates")


# simde library for sim inteercompatibility
include(simde)
target_link_libraries(${PROJECT_NAME} PUBLIC simde)


if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	add_executable(${PROJECT_NAME}_bin test.cpp)
    target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${PROJECT_NAME})
endif()