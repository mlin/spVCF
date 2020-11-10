cmake_minimum_required(VERSION 3.5)
include(ExternalProject)
project(spvcf LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_BUILD_TYPE)

ExternalProject_Add(htslib
    URL https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external
    CONFIGURE_COMMAND bash -c "autoreconf && ./configure --with-libdeflate --disable-libcurl --disable-bz2 --disable-lzma --disable-s3 --disable-gcs"
    PATCH_COMMAND sed -i "s/^CFLAGS .*$/CFLAGS = -O3 -DNDEBUG -march=ivybridge/" Makefile
    BUILD_IN_SOURCE 1
    BUILD_COMMAND bash -c "make -n && make -j$(nproc)"
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    LOG_BUILD ON
  )
ExternalProject_Get_Property(htslib source_dir)
set(HTSLIB_SOURCE_DIR ${source_dir})
ExternalProject_Get_Property(htslib binary_dir)
set(HTSLIB_BINARY_DIR ${binary_dir})

execute_process(COMMAND git describe --tags --long --always
                OUTPUT_VARIABLE GIT_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE)
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGIT_REVISION=\"\\\"${GIT_REVISION}\\\"\"")

add_executable(spvcf src/main.cc src/spVCF.cc src/spVCF.h src/strlcpy.h)
add_dependencies(spvcf htslib)
target_include_directories(spvcf PRIVATE src ${HTSLIB_SOURCE_DIR})
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(spvcf PRIVATE -fdiagnostics-color=auto -march=ivybridge -g)
    set_target_properties(spvcf PROPERTIES LINK_FLAGS "-static-libgcc -static-libstdc++ -pthread")
endif()
target_link_libraries(spvcf ${HTSLIB_BINARY_DIR}/libhts.a libz.a libdeflate.a)

include(CTest)
add_test(NAME tests COMMAND prove -v test/spVCF.t)

# Best practices references:
# https://codingnest.com/basic-cmake/ https://codingnest.com/basic-cmake-part-2/
# https://gist.github.com/mbinna/c61dbb39bca0e4fb7d1f73b0d66a4fd1
