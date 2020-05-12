cmake_minimum_required(VERSION 3.12)

project(obscura)

# Use C++-11
set(CMAKE_CXX_STANDARD 11)
# Require (at least) it
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Don't use e.g. GNU extension (like -std=gnu++11) for portability
set(CMAKE_CXX_EXTENSIONS OFF)

set(BIN_DIR     ${PROJECT_SOURCE_DIR}/bin)
set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(LIB_DIR     ${PROJECT_SOURCE_DIR}/lib)
set(SRC_DIR     ${PROJECT_SOURCE_DIR}/src)
set(TESTS_DIR   ${PROJECT_SOURCE_DIR}/tests)
set(EXTERNAL_DIR  ${PROJECT_SOURCE_DIR}/external)

# IMPORT OF EXTERNAL PROJECTS
include(FetchContent)
# set(FETCHCONTENT_QUIET off)
# libphysica
set(LIBPHYSICA_DIR   ${EXTERNAL_DIR}/libphysica)
FetchContent_Declare(libphysica
  GIT_REPOSITORY  https://github.com/temken/libphysica.git
  GIT_TAG         origin/master
  SOURCE_DIR     "${LIBPHYSICA_DIR}"
  CMAKE_ARGS     "-DCMAKE_BUILD_TYPE=Release -DCODE_COVERAGE=OFF" )
FetchContent_GetProperties(libphysica)
if(NOT libphysica_POPULATED)
    FetchContent_Populate(libphysica)
    add_subdirectory( ${libphysica_SOURCE_DIR} )
endif()

# # libconfig
# set(LIBCONFIG_DIR   ${EXTERNAL_DIR}/libconfig)
# FetchContent_Declare(libconfig
#   GIT_REPOSITORY  https://github.com/hyperrealm/libconfig.git
#   GIT_TAG         origin/master
#   SOURCE_DIR     "${LIBCONFIG_DIR}"
#   CMAKE_ARGS     -DBUILD_EXAMPLES=OFF -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTS=OFF )
# FetchContent_GetProperties(libconfig)
# if(NOT libconfig_POPULATED)
#     FetchContent_Populate(libconfig)
#     add_subdirectory( ${libconfig_SOURCE_DIR} )
# endif()

# Googletest
set(GTEST_DIR   ${EXTERNAL_DIR}/googletest)
FetchContent_Declare(googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.10.0
  SOURCE_DIR     "${GTEST_DIR}/src"
  BINARY_DIR     "${GTEST_DIR}/build"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      "")

# Source and include directories
include_directories( ${INCLUDE_DIR} )
add_subdirectory( ${SRC_DIR} )

# Code Coverage Configuration
if(NOT TARGET coverage_config)
  add_library(coverage_config INTERFACE)
endif()

option(CODE_COVERAGE "Enable coverage reporting" OFF)
if(CODE_COVERAGE)
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O0        # no optimization
    -g         # generate debug info
    --coverage # sets all required flags
  )
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif(CODE_COVERAGE)

option (BUILD_TESTING "Build the testing tree." ON)
# Only build tests if we are the top-level project
# Allows this to be used by super projects with `add_subdirectory`
if (BUILD_TESTING AND (PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR))
  enable_testing()
  FetchContent_GetProperties(googletest)
  if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory(
      ${googletest_SOURCE_DIR} 
      ${googletest_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif()
  add_subdirectory(${TESTS_DIR})
endif()



# #Include my libphys library
# include(ExternalProject)
# set(LIBPHYS_NAME          "libphys-git")
# set(LIBPHYS_GIT           "https://github.com/temken/libphys.git")
# set(LIBPHYS_PREFIX        "${EXPROJ_DIR}/${LIBPHYS_NAME}")
# set(LIBPHYS_SOURCE_DIR    "${EXPROJ_DIR}/${LIBPHYS_NAME}")
# set(LIBPHYS_BUILD_DIR     "${EXPROJ_DIR}/${LIBPHYS_NAME}")
# set(LIBPHYS_INCLUDE_DIR     "${EXPROJ_DIR}/${LIBPHYS_NAME}/include")
# ExternalProject_Add(${LIBPHYS_NAME}
#     PREFIX          ${LIBPHYS_PREFIX}
#     TMP_DIR         ${LIBPHYS_PREFIX}-tmp
#     STAMP_DIR       ${LIBPHYS_PREFIX}-stamp
#     # # - Download Step ------------------
#     GIT_REPOSITORY  ${LIBPHYS_GIT}
#     # - Update Step --------------------
#     UPDATE_COMMAND  ""
#     # - Configure Step -----------------
#     SOURCE_DIR      ${LIBPHYS_SOURCE_DIR}
#     # - Build Step ---------------------
#     BINARY_DIR      ${LIBPHYS_BUILD_DIR}
#     # # - Install Step ------------------
#     # # INSTALL_COMMAND ""
#     CMAKE_ARGS
#     "-DCMAKE_BUILD_TYPE=Release"
# )
# add_library(libphys STATIC IMPORTED)
# set_target_properties(libphys PROPERTIES
#       IMPORTED_LOCATION "${LIBPHYS_SOURCE_DIR}/lib/libphys.a"
#       INTERFACE_INCLUDE_DIRECTORIES "${LIBPHYS_SOURCE_DIR}/src")


# include_directories( ${INCLUDE_DIR} ${LIBPHYS_INCLUDE_DIR})

# add_subdirectory(${SRC_DIR})

# # Code Coverage Configuration
# add_library(coverage_config INTERFACE)
# option(CODE_COVERAGE "Enable coverage reporting" OFF)
# if(CODE_COVERAGE)
#   # Add required flags (GCC & LLVM/Clang)
#   target_compile_options(coverage_config INTERFACE
#     -O0        # no optimization
#     -g         # generate debug info
#     --coverage # sets all required flags
#   )
#   if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
#     target_link_options(coverage_config INTERFACE --coverage)
#   else()
#     target_link_libraries(coverage_config INTERFACE --coverage)
#   endif()
# endif(CODE_COVERAGE)


# option (BUILD_TESTING "Build the testing tree." ON)
# # Only build tests if we are the top-level project
# # Allows this to be used by super projects with `add_subdirectory`
# if (BUILD_TESTING AND (PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR))
#   enable_testing()
#   add_subdirectory(${TESTS_DIR})
# endif()
