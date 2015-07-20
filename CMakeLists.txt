
# === CMake lists for the GaussianLib - (20/07/2015) ===

cmake_minimum_required(VERSION 2.8)
project(GaussianLib)


# === Build path ===

set(OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${OUTPUT_DIR} CACHE PATH "Build directory" FORCE)
set(LIBRARY_OUTPUT_PATH ${OUTPUT_DIR} CACHE PATH "Build directory" FORCE)
set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
set(PROJECT_TEST_DIR "${PROJECT_SOURCE_DIR}/test")


# === Global files ===

file(GLOB_RECURSE FilesAll "${PROJECT_SOURCE_DIR}/include/*.*")


# === Include directories ===

include_directories("${PROJECT_INCLUDE_DIR}")


# === Binaries ===

add_library(gausslib STATIC ${FilesAll})
set_target_properties(gausslib PROPERTIES LINKER_LANGUAGE CXX DEBUG_POSTFIX "D")

add_executable(test1 "${PROJECT_TEST_DIR}/test1.cpp")
set_target_properties(test1 PROPERTIES LINKER_LANGUAGE CXX DEBUG_POSTFIX "D")


