
add_library(loess misc.c loess.c loessf.c loessc.c predict.c supp.c)

include_directories(../f2c)

# Not my code. Supress compiler warnings.
set_target_properties(loess PROPERTIES COMPILE_FLAGS "-w")

