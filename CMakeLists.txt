cmake_minimum_required(VERSION 3.5)
project(StrawmanSparseVCF LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(spvcf src/main.cc src/spVCF.cc src/spVCF.h)
target_include_directories(spvcf PRIVATE src)

include(CTest)
add_test(NAME tests COMMAND prove -v test/spVCF.t)

# Best practices references:
# https://codingnest.com/basic-cmake/ https://codingnest.com/basic-cmake-part-2/
# https://gist.github.com/mbinna/c61dbb39bca0e4fb7d1f73b0d66a4fd1
