
add_subdirectory(gtf)
add_subdirectory(hat-trie)
add_subdirectory(samtools)
add_subdirectory(seqbias)
add_subdirectory(sqlite3)

add_executable(isolator
               common.cpp
               constants.cpp
               emp_dist.cpp
               fragment_model.cpp
               intervals.cpp
               isolator.cpp
               linalg.cpp
               logger.cpp
               pos_table.cpp
               read_set.cpp
               sample_db.cpp
               sampler.cpp
               sparse_mat.cpp
               transcripts.cpp)

include_directories(${isolator_SOURCE_DIR}/src
                    ${Boost_INCLUDE_DIRS}
                    ${ZLIB_INCLUDE_DIRS})

target_link_libraries(isolator
                      gtf
                      hat-trie
                      samtools
                      seqbias
                      sqlite3
                      ${Boost_SYSTEM_LIBRARY}
                      ${Boost_THREAD_LIBRARY}
                      ${Boost_TIMER_LIBRARY}
                      ${ZLIB_LIBRARIES}
                      ${GSL_LIBRARIES})

