
add_subdirectory(cblas)
add_subdirectory(f2c)
add_subdirectory(gtf)
add_subdirectory(hat-trie)
add_subdirectory(loess)
add_subdirectory(samtools)
add_subdirectory(seqbias)
add_subdirectory(sqlite3)
add_subdirectory(stan)

add_executable(isolator
               analysis.cpp
               analyze.cpp
               common.cpp
               constants.cpp
               cpuid.cpp
               dirichlet.c
               emp_dist.cpp
               fragment_model.cpp
               gc.cpp
               intervals.cpp
               isolator.cpp
               linalg.cpp
               logger.cpp
               pos_table.cpp
               read_set.cpp
               sample_db.cpp
               sampler.cpp
               sparse_mat.cpp
               switchtest.cpp
               transcripts.cpp)

include_directories(${isolator_SOURCE_DIR}/src
                    ${Boost_INCLUDE_DIRS}
                    ${ZLIB_INCLUDE_DIRS})

target_link_libraries(isolator
                      loess
                      cblas
                      f2c
                      gtf
                      hat-trie
                      samtools
                      seqbias
                      sqlite3
                      stan
                      nlopt
                      ${Boost_SYSTEM_LIBRARY}
                      ${Boost_THREAD_LIBRARY}
                      ${Boost_TIMER_LIBRARY}
                      ${ZLIB_LIBRARIES}
                      ${GSL_LIBRARIES}
                      ${CMAKE_THREAD_LIBS_INIT})

add_custom_command(
    OUTPUT "${CMAKE_SOURCE_DIR}/analyze_model.hpp"
    COMMAND stanc --no_main --o=${CMAKE_SOURCE_DIR}/src/analyze_model.hpp analyze.stan
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/src"
    DEPENDS analyze.stan)

#set_source_files_properties(analyze.cpp
                            #PROPERTIES OBJECT_DEPENDS analyze_model.cpp)

include_directories("eigen")
