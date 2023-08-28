cmake_minimum_required(VERSION 3.10)

# set the project name
project(indirectPredicates)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add the executable
add_executable(${PROJECT_NAME} 
    test.cpp
)

# Compiler-specific options
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME} PUBLIC "/fp:strict")
	# use intrinsic functions
	target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")
	# reserve enough stack size
	target_link_options(${PROJECT_NAME} PUBLIC "/STACK:8421376")
	# turn off annoying warnings
	target_compile_options(${PROJECT_NAME} PUBLIC "/D _CRT_SECURE_NO_WARNINGS")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	# set standard optimization level
	target_compile_options(${PROJECT_NAME} PUBLIC -O2)
	# reserve enough stack size	
	target_compile_options(${PROJECT_NAME} PUBLIC -Wl,-z,stacksize=8421376)
	# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME} PUBLIC -frounding-math)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	# reserve enough stack size	
	target_compile_options(${PROJECT_NAME} PUBLIC -Wl,-z,stacksize=8421376)
	# Disable all optimizations
	target_compile_options(${PROJECT_NAME} PUBLIC -O0)
endif()

# Public include directory
target_include_directories(${PROJECT_NAME} PUBLIC 
	.
	include
)
