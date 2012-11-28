
add_library(sqlite3 sqlite3.c)

add_definitions(-DSQLITE_OMIT_LOAD_EXTENSION)
