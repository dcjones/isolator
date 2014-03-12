
add_library(yaml
            api.c
            dumper.c
            emitter.c
            loader.c
            parser.c
            reader.c
            scanner.c
            writer.c)

# Not my code. Supress compiler warnings.
set_target_properties(yaml PROPERTIES COMPILE_FLAGS "-w")

include_directories()
