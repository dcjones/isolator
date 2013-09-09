cmake_minimum_required(VERSION 2.6)
project(isolator)

include(CheckIncludeFile)

set(CMAKE_COLOR_MAKEFILE ON)

execute_process(COMMAND sh ${CMAKE_SOURCE_DIR}/git-version-gen
                OUTPUT_VARIABLE VERSION
                OUTPUT_STRIP_TRAILING_WHITESPACE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

set(Boost_USE_MULTITHREADED ON)
find_package(Boost 1.48.0 COMPONENTS thread system timer REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Threads)

# Check for SSE or AVX intrinsics
check_include_file(immintrin.h HAVE_IMMINTRIN_H)

set(ARCH_FLAGS "-msse2 -msse4.1 -msse4.2 -mavx")
set(CMAKE_C_FLAGS "-std=gnu99 -g -O3 -Wall -Wextra -DNDEBUG ${ARCH_FLAGS}")
set(CMAKE_CXX_FLAGS "-std=gnu++98 -g -O3 -Wall -Wextra -DNDEBUG ${ARCH_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG "-std=gnu99 -g -O0 -Wall -Wextra ${ARCH_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "-std=gnu++98 -g -O0 -Wall -Wextra ${ARCH_FLAGS}")


configure_file("${CMAKE_SOURCE_DIR}/src/config.h.in"
               "${CMAKE_SOURCE_DIR}/src/config.h")

add_subdirectory(src)

