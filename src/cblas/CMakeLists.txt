
add_library(cblas
    ddot.c daxpy.c dnrm2.c dqrsl.c dcopy.c dscal.c dswap.c dsvdc.c drot.c
    drotg.c)

# Not my code. Supress compiler warnings.
set_target_properties(cblas PROPERTIES COMPILE_FLAGS "-w")

include_directories(../f2c)

