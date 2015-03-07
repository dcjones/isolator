
include_directories(${isolator_SOURCE_DIR}/src
                    ${Boost_INCLUDE_DIRS}
                    ${HDF5_INCLUDE_DIRS}
                    ${ZLIB_INCLUDE_DIRS}
                    ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(gtf)
add_subdirectory(hat-trie)
add_subdirectory(nlopt)
add_subdirectory(samtools)
add_subdirectory(seqbias)
add_subdirectory(yaml)

add_executable(isolator
               analyze.cpp
               common.cpp
               constants.cpp
               cpuid.cpp
               emp_dist.cpp
               fastmath.cpp
               fastmath-vanilla.cpp
               fastmath-sse2.cpp
               fastmath-sse4.cpp
               fastmath-avx.cpp
               fragbias.cpp
               fragment_model.cpp
               gcbias.cpp
               hdf5.cpp
               intervals.cpp
               isolator.cpp
               logger.cpp
               pos_table.cpp
               read_set.cpp
               sampler.cpp
               shredder.cpp
               tpbias.cpp
               summarize.cpp
               transcripts.cpp)

SET_SOURCE_FILES_PROPERTIES(fastmath-avx.cpp PROPERTIES COMPILE_FLAGS "-msse2 -msse4.1 -msse4.2 -mavx")
SET_SOURCE_FILES_PROPERTIES(fastmath-sse4.cpp PROPERTIES COMPILE_FLAGS "-msse2 -msse4.1 -msse4.2")
SET_SOURCE_FILES_PROPERTIES(fastmath-sse2.cpp PROPERTIES COMPILE_FLAGS "-msse2")

target_link_libraries(isolator
                      gtf
                      hat-trie
                      samtools
                      seqbias
                      yaml
                      nlopt
                      ${Boost_SYSTEM_LIBRARY}
                      ${Boost_THREAD_LIBRARY}
                      ${Boost_TIMER_LIBRARY}
                      ${Boost_CHRONO_LIBRARY}
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
                      nlopt
                      ${Boost_SYSTEM_LIBRARY}
                      ${Boost_THREAD_LIBRARY}
                      ${Boost_TIMER_LIBRARY}
                      ${Boost_CHRONO_LIBRARY}
                      ${ZLIB_LIBRARIES}
                      ${CMAKE_THREAD_LIBS_INIT})

if(UNIX AND NOT APPLE)
  target_link_libraries(isolator rt)
  target_link_libraries(samcnt rt)
endif()


install(TARGETS isolator RUNTIME DESTINATION bin)

# Fancy automatic git versions
add_custom_target(gitversion_target
                  ALL DEPENDS gitversion)

add_custom_command(OUTPUT gitversion
                   COMMAND ${CMAKE_COMMAND}
                   -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
                   -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/gitversion.cmake)

add_dependencies(isolator gitversion_target)


