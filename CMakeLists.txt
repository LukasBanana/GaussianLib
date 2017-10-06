
# === CMake lists for the GaussianLib - (20/07/2015) ===

cmake_minimum_required(VERSION 2.8)
project(GaussianLib)


# === Build path ===

set(OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${OUTPUT_DIR} CACHE PATH "Build directory" FORCE)
set(LIBRARY_OUTPUT_PATH ${OUTPUT_DIR} CACHE PATH "Build directory" FORCE)
set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
set(PROJECT_TEST_DIR "${PROJECT_SOURCE_DIR}/test")


# === Options ===

option(GAUSSLIB_ENABLE_SSE "Enable SSE support" OFF)
option(GAUSSLIB_ENABLE_SSE2 "Enable SSE2 support" OFF)

if(GAUSSLIB_ENABLE_SSE)
	add_definitions(-DGS_ENABLE_SSE)
endif()

if(GAUSSLIB_ENABLE_SSE2)
	add_definitions(-DGS_ENABLE_SSE2)
endif()

if(UNIX AND (GAUSSLIB_ENABLE_SSE OR GAUSSLIB_ENABLE_SSE2))
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
endif()


# === Global files ===

file(GLOB_RECURSE FilesAll "${PROJECT_SOURCE_DIR}/include/*.*")


# === Include directories ===

include_directories("${PROJECT_INCLUDE_DIR}")


# === Binaries ===

add_library(gausslib STATIC ${FilesAll})
set_target_properties(gausslib PROPERTIES LINKER_LANGUAGE CXX DEBUG_POSTFIX "D")
target_compile_features(gausslib PRIVATE cxx_range_for)

add_executable(test1 "${PROJECT_TEST_DIR}/test1.cpp")
set_target_properties(test1 PROPERTIES LINKER_LANGUAGE CXX DEBUG_POSTFIX "D")
target_compile_features(test1 PRIVATE cxx_range_for)


