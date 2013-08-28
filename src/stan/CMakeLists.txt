
add_library(stan
            agrad/rev/var_stack.cpp
            math/matrix.cpp
            agrad/matrix.cpp)

include_directories(".."
                    "../eigen")

set_target_properties(stan PROPERTIES COMPILE_FLAGS "-w")
