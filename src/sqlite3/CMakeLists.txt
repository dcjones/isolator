
add_library(sqlite3 sqlite3.c)

add_definitions(-DSQLITE_OMIT_LOAD_EXTENSION)

# Not my code. Supress compiler warnings.
set_target_properties(sqlite3 PROPERTIES COMPILE_FLAGS "-w")

