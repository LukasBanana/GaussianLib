
# === CMake lists for the GaussianLib - (20/07/2015) ===

cmake_minimum_required(VERSION 3.19)
project(GaussianLib)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# === Build path ===

set(OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${OUTPUT_DIR} CACHE PATH "Build directory" FORCE)
set(LIBRARY_OUTPUT_PATH ${OUTPUT_DIR} CACHE PATH "Build directory" FORCE)
set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
set(PROJECT_TEST_DIR "${PROJECT_SOURCE_DIR}/test")


# === Options ===

option(GaussLib_REAL_DOUBLE "Use double precision floating-points for real types" OFF)
option(GaussLib_ASSERT_EXCEPTION "Enable runtime exceptions for assertions" OFF)
option(GaussLib_ENABLE_SWIZZLE_OPERATOR "Enable swizzle operator for Vector class" OFF)
option(GaussLib_ENABLE_INVERSE_OPERATOR "Enable inverse matrix operator (A ^ -1)" OFF)
option(GaussLib_DISABLE_AUTO_INIT "Disable automatic initialization" OFF)
option(GaussLib_ROW_MAJOR_STORAGE "Use row-major storage (column-major storage otherwise)" OFF)
option(GaussLib_ROW_VECTORS "Use row-vectors (column-vectors otherwise)" OFF)


# === Macros ===

if(GaussLib_REAL_DOUBLE)
	add_definitions(-DGS_REAL_DOUBLE)
endif()

if(GaussLib_ASSERT_EXCEPTION)
	add_definitions(-DGS_ASSERT_EXCEPTION)
endif()

if(GaussLib_ENABLE_SWIZZLE_OPERATOR)
	add_definitions(-DGS_ENABLE_SWIZZLE_OPERATOR)
endif()

if(GaussLib_ENABLE_INVERSE_OPERATOR)
	add_definitions(-DGS_ENABLE_INVERSE_OPERATOR)
endif()

if(GaussLib_DISABLE_AUTO_INIT)
	add_definitions(-DGS_DISABLE_AUTO_INIT)
endif()

if(GaussLib_ROW_MAJOR_STORAGE)
	add_definitions(-DGS_ROW_MAJOR_STORAGE)
endif()

if(GaussLib_ROW_VECTORS)
	add_definitions(-DGS_ROW_VECTORS)
endif()


# === Global files ===

file(GLOB_RECURSE FilesAll "${PROJECT_SOURCE_DIR}/include/*.*")


# === Include directories ===

#include_directories(${PROJECT_INCLUDE_DIR})


# === Binaries ===

add_library(gausslib INTERFACE ${FilesAll})
target_include_directories(gausslib INTERFACE ${PROJECT_INCLUDE_DIR})

set(
	FilesTest1
	${PROJECT_TEST_DIR}/test1.h
	${PROJECT_TEST_DIR}/test1.cpp
	${PROJECT_TEST_DIR}/test1_main.cpp
)

add_executable(test1 ${FilesTest1})
target_link_libraries(test1 PRIVATE gausslib)

