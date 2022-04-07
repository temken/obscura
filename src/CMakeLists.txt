# Executable
add_executable(obscura
    main.cpp )
     
if(BUILD_SHARED_LIBS)
    set_target_properties(obscura PROPERTIES
        INSTALL_RPATH "${LIB_DIR}") # required if libobscura is built as shared library
endif()

target_compile_options(obscura PUBLIC -Wall -pedantic)

target_link_libraries(obscura
    PUBLIC 
        coverage_config
        libobscura )

target_include_directories(obscura
    PRIVATE
        ${GENERATED_DIR} )

install(TARGETS obscura DESTINATION ${BIN_DIR})

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

# Library (default is static, but can be build as shared using -DBUILD_SHARED_LIBS=ON in the configuration step)
add_library(libobscura
    Astronomy.cpp
    Configuration.cpp
    Direct_Detection.cpp
    Direct_Detection_ER.cpp
    Direct_Detection_Ionization.cpp
    Direct_Detection_Migdal.cpp
    Direct_Detection_Nucleus.cpp
    Direct_Detection_Crystal.cpp
    DM_Distribution.cpp
    DM_Halo_Models.cpp
    DM_Particle.cpp
    DM_Particle_Standard.cpp
    Experiments.cpp
    Target_Atom.cpp
    Target_Crystal.cpp
    Target_Nucleus.cpp
)

set_target_properties(libobscura PROPERTIES PREFIX "")
if(BUILD_SHARED_LIBS)
    set_target_properties(libobscura PROPERTIES SUFFIX ".so")
endif()

target_compile_options(libobscura PUBLIC -Wall -pedantic)

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
