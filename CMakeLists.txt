cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(FemGauntlet LANGUAGES C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(BUILD_SHARED_LIBS ON)

include_directories(.)

# Add the directory with the mesh generator
add_subdirectory(genmesh)

# Add the directory with your submission
add_subdirectory(project)


add_executable(
    myFem
    ./benchmark.c
)

find_library(m m)
target_link_libraries(
    myFem 
    PUBLIC femsolver 
)
target_link_libraries(myFem PUBLIC m)

# copy libfemsolve.dll to the build directory if on windows
if(WIN32)
    add_custom_command(
        TARGET myFem POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
            $<TARGET_FILE:femsolver>
            $<TARGET_FILE_DIR:myFem>
    )
endif()
