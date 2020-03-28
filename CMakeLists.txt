include_guard()

cmake_minimum_required (VERSION 3.5.1)

project (Linear CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/release)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_BINARY_DIR}/releasewd)

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
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
		target_compile_options(Linear PUBLIC
			-Weverything
			-Wno-c++98-compat
			-Wno-c++98-compat-pedantic)
	elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		target_compile_options(Linear PUBLIC -Wall)
	endif()
else()
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
		target_compile_options(Linear PUBLIC -W -Wall -o3)
	elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		target_compile_options(Linear PUBLIC -Wall -o3)
	endif()
ENDIF()

include_directories(./includes/)

enable_testing ()
add_subdirectory (tests)