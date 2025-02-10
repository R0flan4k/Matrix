cmake_minimum_required(VERSION 3.11)
project(Matrix)
if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        set(CMAKE_CXX_FLAGS_DEBUG "-fsanitize=address -fsanitize=undefined")
endif()

set(SRCS
	${CMAKE_SOURCE_DIR}/Matrix/src/matrix.cpp
)
add_executable(Matrix.x)

target_compile_features(Matrix.x PUBLIC cxx_std_20)
if((NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC") AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
        target_compile_options(Matrix.x PUBLIC -std=c++20 -Wall -g -O0)
endif()
target_sources(Matrix.x PRIVATE ${SRCS})
target_include_directories(Matrix.x PUBLIC "${CMAKE_SOURCE_DIR}/Matrix/include")

set(CLANG_FORMAT_SRCS
        ${SRCS}
	${CMAKE_SOURCE_DIR}/Matrix/include/*.h
        ${CMAKE_SOURCE_DIR}/unit_tests/*.cpp
)
add_custom_target(
        clangformat
        ALL
        COMMAND clang-format -style=file -i ${CLANG_FORMAT_SRCS}
)

add_subdirectory(${CMAKE_SOURCE_DIR}/unit_tests)
add_subdirectory(${CMAKE_SOURCE_DIR}/e2e)

