
add_library(cblas
    ddot.c daxpy.c dnrm2.c dqrsl.c dcopy.c dscal.c dswap.c dsvdc.c drot.c
    drotg.c)

include_directories(../f2c)

