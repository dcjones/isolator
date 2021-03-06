cmake_minimum_required(VERSION 2.6)
project(isolator)

include(CheckIncludeFile)
include(CheckFunctionExists)

set(CMAKE_COLOR_MAKEFILE ON)

#execute_process(COMMAND sh ${CMAKE_SOURCE_DIR}/git-version-gen
                #OUTPUT_VARIABLE VERSION
                #OUTPUT_STRIP_TRAILING_WHITESPACE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

set(Boost_USE_MULTITHREADED ON)
find_package(Boost 1.48.0 COMPONENTS thread system timer chrono REQUIRED)
find_package(HDF5 COMPONENTS C HL REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Threads)

check_function_exists(isnan HAVE_ISNAN)
check_include_file(stdint.h HAVE_STDINT_H)

# Check for SSE or AVX intrinsics
check_include_file(immintrin.h HAVE_IMMINTRIN_H)

if(BUILD_STATIC_EXECUTABLES)
    set(CMAKE_EXE_LINKER_FLAGS -static)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
endif(BUILD_STATIC_EXECUTABLES)

set(CMAKE_C_FLAGS_RELEASE "-std=gnu99 -g -O3 -Wall -Wextra -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-std=gnu++98 -g -O3 -Wall -Wextra -DNDEBUG")
set(CMAKE_C_FLAGS_DEBUG "-std=gnu99 -g -O0 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-std=gnu++98 -g -O0 -Wall -Wextra")


configure_file("${CMAKE_SOURCE_DIR}/src/config.h.in"
               "${CMAKE_SOURCE_DIR}/src/config.h")

add_subdirectory(src)


