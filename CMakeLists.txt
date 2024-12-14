cmake_minimum_required(VERSION 3.11)
project(Matrix)
set(CMAKE_CXX_FLAGS "-O2 -std=c++20")
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g -O0 -fsanitize=address -fsanitize=undefined")

set(SRCS
	${CMAKE_SOURCE_DIR}/Matrix/src/matrix.cpp
)
add_executable(Matrix.x)

target_sources(Matrix.x PRIVATE ${SRCS})
target_include_directories(Matrix.x PUBLIC "${CMAKE_SOURCE_DIR}/Matrix/include")

set(CLANG_FORMAT_SRCS
        ${SRCS}
	${CMAKE_SOURCE_DIR}/Matrix/include/*.h
        ${CMAKE_SOURCE_DIR}/tests/*.cpp
)
add_custom_target(
        clangformat
        ALL
        COMMAND clang-format -style=file -i ${CLANG_FORMAT_SRCS}
)


add_subdirectory(${CMAKE_SOURCE_DIR}/tests)

