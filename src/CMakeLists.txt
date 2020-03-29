add_library(darkmatter STATIC
    Astronomy.cpp
    DM_Distribution.cpp
    DM_Particle.cpp
    Target_Nucleus.cpp
)

target_include_directories(darkmatter  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_options(darkmatter PUBLIC -Wall -pedantic)

# Include code-coverage settings:
target_link_libraries(darkmatter PUBLIC coverage_config)

install(TARGETS darkmatter DESTINATION ${LIB_DIR})


# file(GLOB SOURCES ${SRC_DIR}/*.cpp)

# add_executable(DM_Detection ${SOURCES})

add_executable(DM_Detection main.cpp)

target_link_libraries(DM_Detection
	PUBLIC darkmatter
	PUBLIC coverage_config 
	PUBLIC libphys)

# add_test(NAME DM_Detection_test COMMAND DM_Detection)
# add_test(NAME DM_Detection_test COMMAND DM_Detection WORKING_DIRECTORY ${BIN_DIR})

install(TARGETS DM_Detection DESTINATION ${BIN_DIR})


