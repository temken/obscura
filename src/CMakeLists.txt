add_library(obscura STATIC
    Astronomy.cpp
    Configuration.cpp
    Direct_Detection.cpp
    Direct_Detection_Ionization.cpp
    Direct_Detection_Nucleus.cpp
    Direct_Detection_Semiconductor.cpp
    DM_Distribution.cpp
    DM_Particle.cpp
    DM_Particle_Standard.cpp
    Target_Electron.cpp
    Target_Nucleus.cpp
)

target_compile_options(obscura PUBLIC -Wall -pedantic)

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

target_include_directories(obscura 
                        PUBLIC
                            ${LIBCONFIG_INCLUDE_DIRs})

# Link libraries:
target_link_libraries(obscura 
                        PUBLIC
                            coverage_config 
                            ${LIBCONFIGPP_LIBRARY})

install(TARGETS obscura DESTINATION ${LIB_DIR})


add_executable(obscura-bin main.cpp)
SET_TARGET_PROPERTIES(obscura-bin
  PROPERTIES OUTPUT_NAME obscura)

target_link_libraries(obscura-bin
	PUBLIC obscura
	PUBLIC coverage_config 
	PUBLIC libphys)

install(TARGETS obscura-bin DESTINATION ${BIN_DIR})


