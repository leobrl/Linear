include_guard()

set( CMAKE_VERBOSE_MAKEFILE on )

cmake_minimum_required (VERSION 3.5.1)

project (Linear CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/release)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_BINARY_DIR}/release_wd)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL ${CMAKE_BINARY_DIR}/min_size_rel)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED 17)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	set(CMAKE_CXX_FLAGS_DEBUG "-ggdb3 -O0" CACHE STRING "" FORCE)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "-ggdb3 -O0" CACHE STRING "" FORCE)
else()
    message(FATAL " Compiler not supported.")
endif()

add_executable(Linear src/Linear_main.cpp)

IF(CMAKE_BUILD_TYPE MATCHES "Debug")
	message(STATUS "Target compiled with debug compiler options.")
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
		target_compile_options(Linear PUBLIC
			-Weverything
			-Wno-c++98-compat
			-Wno-c++98-compat-pedantic)
	elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		target_compile_options(Linear PUBLIC -Wall)
	endif()
else()
	message(STATUS "Target compiled with release compiler options.")
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
		target_compile_options(Linear PUBLIC -W -Wall -O3)
	elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		target_compile_options(Linear PUBLIC -Wall -O3)
	endif()
ENDIF()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(Linear PUBLIC OpenMP::OpenMP_CXX)
else()
	message(FATAL "OpenMP not found.")
endif()

include_directories(./includes/)

enable_testing ()
add_subdirectory (tests)