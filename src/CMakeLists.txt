
include_directories(${isolator_SOURCE_DIR}/src
                    ${Boost_INCLUDE_DIRS}
                    ${HDF5_INCLUDE_DIRS}
                    ${ZLIB_INCLUDE_DIRS})

add_subdirectory(cblas)
add_subdirectory(f2c)
add_subdirectory(gtf)
add_subdirectory(hat-trie)
add_subdirectory(loess)
add_subdirectory(samtools)
add_subdirectory(seqbias)

add_executable(isolator
               analyze.cpp
               common.cpp
               constants.cpp
               cpuid.cpp
               emp_dist.cpp
               fragment_model.cpp
               gc.cpp
               intervals.cpp
               isolator.cpp
               linalg.cpp
               logger.cpp
               pos_table.cpp
               read_set.cpp
               sampler.cpp
               shredder.cpp
               summarize.cpp
               sparse_mat.cpp
               transcripts.cpp)

target_link_libraries(isolator
                      loess
                      cblas
                      f2c
                      gtf
                      hat-trie
                      samtools
                      seqbias
                      ${Boost_SYSTEM_LIBRARY}
                      ${Boost_THREAD_LIBRARY}
                      ${Boost_TIMER_LIBRARY}
                      ${HDF5_LIBRARIES}
                      ${ZLIB_LIBRARIES}
                      ${CMAKE_THREAD_LIBS_INIT})

