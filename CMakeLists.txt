cmake_minimum_required(VERSION 3.2)

project(DM_Detection)

# Compiler settings
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Directories
set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(BIN_DIR ${PROJECT_SOURCE_DIR}/bin)
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(EXPROJ_DIR    "${CMAKE_HOME_DIRECTORY}/external")
set(TESTS_DIR ${PROJECT_SOURCE_DIR}/tests)
set(LIB_DIR ${PROJECT_SOURCE_DIR}/lib)

#Include my libphys library
include(ExternalProject)
set(LIBPHYS_NAME          "libphys-git")
set(LIBPHYS_GIT           "https://github.com/temken/libphys.git")
set(LIBPHYS_PREFIX        "${EXPROJ_DIR}/${LIBPHYS_NAME}")
set(LIBPHYS_SOURCE_DIR    "${EXPROJ_DIR}/${LIBPHYS_NAME}")
set(LIBPHYS_BUILD_DIR     "${EXPROJ_DIR}/${LIBPHYS_NAME}")
set(LIBPHYS_INCLUDE_DIR     "${EXPROJ_DIR}/${LIBPHYS_NAME}/include")
ExternalProject_Add(${LIBPHYS_NAME}
    PREFIX          ${LIBPHYS_PREFIX}
    TMP_DIR         ${LIBPHYS_PREFIX}-tmp
    STAMP_DIR       ${LIBPHYS_PREFIX}-stamp
    # # - Download Step ------------------
    GIT_REPOSITORY  ${LIBPHYS_GIT}
    # - Update Step --------------------
    UPDATE_COMMAND  ""
    # - Configure Step -----------------
    SOURCE_DIR      ${LIBPHYS_SOURCE_DIR}
    # - Build Step ---------------------
    BINARY_DIR      ${LIBPHYS_BUILD_DIR}
    # # - Install Step ------------------
    # # INSTALL_COMMAND ""
    CMAKE_ARGS
    "-DCMAKE_BUILD_TYPE=Release"
)
add_library(libphys STATIC IMPORTED)
set_target_properties(libphys PROPERTIES
      IMPORTED_LOCATION "${LIBPHYS_SOURCE_DIR}/lib/libphys.a"
      INTERFACE_INCLUDE_DIRECTORIES "${LIBPHYS_SOURCE_DIR}/src")


include_directories( ${INCLUDE_DIR} ${LIBPHYS_INCLUDE_DIR})

add_subdirectory(${SRC_DIR})

# Code Coverage Configuration
add_library(coverage_config INTERFACE)
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
  add_subdirectory(${TESTS_DIR})
endif()
