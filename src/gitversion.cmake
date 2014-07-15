
execute_process(COMMAND sh ${PROJECT_SOURCE_DIR}/git-version-gen
                OUTPUT_VARIABLE GITVERSION
                OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND sh ${PROJECT_SOURCE_DIR}/git-commit-gen
                OUTPUT_VARIABLE GITCOMMIT
                OUTPUT_STRIP_TRAILING_WHITESPACE)

file(WRITE _gitversion.hpp
     "#ifndef ISOLATOR_GITVERSION_HPP\n#define ISOLATOR_GITVERSION_HPP\n#define GITVERSION \"${GITVERSION}\"\n#define GITCOMMIT \"${GITCOMMIT}\"\n#endif\n")

execute_process(COMMAND ${CMAKE_COMMAND}
                -E copy_if_different
                _gitversion.hpp
                gitversion.hpp)

