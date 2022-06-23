# Set the compiler as described here: https://cmake.org/Wiki/CMake_FAQ#How_do_I_use_a_different_compiler.3F
#
#  Verify that Cmake is up to date
#
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
#SET(CMAKE_LEGACY_CYGWIN_WIN32 0)

#
#  Project name
#
PROJECT("pysces")

#
#  User-defined Find modules
#
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

#
#  Cmake configuration
#
MESSAGE("Setting up the CMake configuration...")
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "RELEASE")
ENDIF()

#
# Now building the project
#

#
#  Copy python files
#
add_custom_target(pysces ALL
    COMMAND cmake -E remove_directory ${CMAKE_BINARY_DIR}/src
    COMMAND cmake -E copy_directory ${CMAKE_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/src
    COMMAND cmake -E copy ${CMAKE_SOURCE_DIR}/src/__init__.py ${CMAKE_BINARY_DIR}/src/__init__.py
)





