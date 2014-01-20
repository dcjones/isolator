
include_directories(${isolator_SOURCE_DIR}/src
                    ${Boost_INCLUDE_DIRS}
                    ${HDF5_INCLUDE_DIRS}
                    ${ZLIB_INCLUDE_DIRS}
                    ${CMAKE_CURRENT_BINARY_DIR})

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
               hdf5.cpp
               intervals.cpp
               isolator.cpp
               linalg.cpp
               logger.cpp
               pos_table.cpp
               read_set.cpp
               sampler.cpp
               shredder.cpp
               summarize.cpp
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

add_executable(samcnt
               samcnt.cpp
               common.cpp
               constants.cpp
               intervals.cpp
               logger.cpp
               transcripts.cpp
               read_set.cpp)

target_link_libraries(samcnt
                      hat-trie
                      samtools
                      seqbias
                      gtf
                      ${Boost_SYSTEM_LIBRARY}
                      ${Boost_THREAD_LIBRARY}
                      ${Boost_TIMER_LIBRARY}
                      ${ZLIB_LIBRARIES}
                      ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS isolator RUNTIME DESTINATION bin)

# Fancy automatic git versions
add_custom_target(gitversion_target
                  ALL DEPENDS gitversion)

add_custom_command(OUTPUT gitversion
                   COMMAND ${CMAKE_COMMAND}
                   -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
                   -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/gitversion.cmake)

set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/gitversion.hpp
                            PROPERTIES GENERATED TRUE
                            HEADER_FILE_ONLY TRUE)

add_dependencies(isolator gitversion_target)


