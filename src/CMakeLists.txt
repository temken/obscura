# Executable
add_executable(obscura
    main.cpp)

if(BUILD_SHARED_LIBS)
    set_target_properties(obscura PROPERTIES
        INSTALL_RPATH "${LIB_DIR}") # required if libobscura is built as shared library
endif()

target_compile_options(obscura PUBLIC -Wall -pedantic)

target_link_libraries(obscura
    PUBLIC
    coverage_config
    libobscura)

target_include_directories(obscura
    PRIVATE
    ${GENERATED_DIR})

install(TARGETS obscura DESTINATION ${BIN_DIR})

# Library (default is static, but can be build as shared using -DBUILD_SHARED_LIBS=ON in the configuration step)
file(GLOB FILES "*.cpp")
list(REMOVE_ITEM FILES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_library(libobscura
    ${FILES})

set_target_properties(libobscura PROPERTIES PREFIX "")

if(BUILD_SHARED_LIBS)
    set_target_properties(libobscura PROPERTIES SUFFIX ".so")
endif()

target_compile_options(libobscura PUBLIC -Wall -pedantic)

target_include_directories(libobscura
    PRIVATE
    ${GENERATED_DIR}
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${INCLUDE_DIR})

target_link_libraries(libobscura
    PUBLIC
    coverage_config
    libphysica
)

install(TARGETS libobscura DESTINATION ${LIB_DIR})
