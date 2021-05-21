# Executable
add_executable(obscura
    main.cpp )

target_compile_options(obscura PUBLIC -Wall -pedantic)

target_link_libraries(obscura
    PUBLIC 
        coverage_config
        libobscura )

target_include_directories(obscura
    PRIVATE
        ${GENERATED_DIR} )

install(TARGETS obscura DESTINATION ${BIN_DIR})

# Static library
add_library(libobscura STATIC
    Astronomy.cpp
    Configuration.cpp
    Direct_Detection.cpp
    Direct_Detection_ER.cpp
    Direct_Detection_Ionization.cpp
    Direct_Detection_Migdal.cpp
    Direct_Detection_Nucleus.cpp
    Direct_Detection_Crystal.cpp
    DM_Distribution.cpp
    DM_Particle.cpp
    DM_Particle_Standard.cpp
    Experiments.cpp
    Target_Atom.cpp
    Target_Crystal.cpp
    Target_Nucleus.cpp
)

set_target_properties(libobscura PROPERTIES PREFIX "")

target_compile_options(libobscura PUBLIC -Wall -pedantic)

# Find and include libconfig
find_path(LIBCONFIG_INCLUDE_DIRs libconfig.h++
    /usr/local/include
    /usr/include
    ~/libs/include
    ~/include )
find_library(LIBCONFIGPP_LIBRARY config++
    /usr/local/lib
    /usr/lib
    ~/libs/lib
    ~/lib )

target_include_directories(libobscura 
    PRIVATE
        ${GENERATED_DIR}
     PUBLIC
        ${LIBCONFIG_INCLUDE_DIRs}
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${INCLUDE_DIR} )

target_link_libraries(libobscura
    PUBLIC
        coverage_config 
        libphysica
        ${LIBCONFIGPP_LIBRARY} )

install(TARGETS libobscura DESTINATION ${LIB_DIR})

