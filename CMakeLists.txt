cmake_minimum_required(VERSION 3.5)
project(StrawmanSparseVCF LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

execute_process(COMMAND git describe --tags --long --dirty --always
                OUTPUT_VARIABLE GIT_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE)
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGIT_REVISION=\"\\\"${GIT_REVISION}\\\"\"")

add_executable(spvcf src/main.cc src/spVCF.cc src/spVCF.h)
target_include_directories(spvcf PRIVATE src)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(spvcf PRIVATE -fdiagnostics-color=auto -march=ivybridge -g)
    set_target_properties(spvcf PROPERTIES LINK_FLAGS -static-libstdc++)
endif()

include(CTest)
add_test(NAME tests COMMAND prove -v test/spVCF.t)

# Best practices references:
# https://codingnest.com/basic-cmake/ https://codingnest.com/basic-cmake-part-2/
# https://gist.github.com/mbinna/c61dbb39bca0e4fb7d1f73b0d66a4fd1
