# empty variable to keep track of the `make example` dependencies
set(examples_depends)

add_subdirectory(programs)

# This is the only one requiring ROOT and BAT

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT)
if(ROOT_FOUND)

    find_package(BAT)
    if(BAT_FOUND)
        add_subdirectory(bat_gen)
    else(BAT_FOUND)
        message(WARNING "BAT not found: bat_gen examples won't be compilable.")
    endif(BAT_FOUND)

else(ROOT_FOUND)
    message(WARNING "ROOT not found: bat_gen examples won't be compilable.")
endif(ROOT_FOUND)

# Create custom 'examples' target
add_custom_target(examples DEPENDS ${examples_depends})
