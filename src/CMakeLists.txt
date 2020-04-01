add_library(darkmatter STATIC
    Astronomy.cpp
    Configuration.cpp
    Direct_Detection.cpp
    Direct_Detection_Nucleus.cpp
    DM_Distribution.cpp
    DM_Particle.cpp
    DM_Particle_Standard.cpp
    Target_Nucleus.cpp
)

target_compile_options(darkmatter PUBLIC -Wall -pedantic)

# Find and include libconfig
find_path(LIBCONFIG_INCLUDE_DIRs libconfig.h++
/usr/local/include
/usr/include
~/libs/include
~/include
)
find_library(LIBCONFIGPP_LIBRARY config++
/usr/local/lib
/usr/lib
~/libs/lib
)

target_include_directories(darkmatter 
                        PUBLIC
                            ${LIBCONFIG_INCLUDE_DIRs})

# Link libraries:
target_link_libraries(darkmatter 
                        PUBLIC
                            coverage_config 
                            ${LIBCONFIGPP_LIBRARY})

install(TARGETS darkmatter DESTINATION ${LIB_DIR})


add_executable(DM_Detection main.cpp)

target_link_libraries(DM_Detection
	PUBLIC darkmatter
	PUBLIC coverage_config 
	PUBLIC libphys)

install(TARGETS DM_Detection DESTINATION ${BIN_DIR})


