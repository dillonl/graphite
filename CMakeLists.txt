# =================================
# gwiz
#
# gwiz/CMakeLists.tx
# =================================

#Set up min CMake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


#Set project name
PROJECT(gwiz)

SET(VERSION_MAJOR "0")
SET(VERSION_MINOR "1")
SET(VERSION_PATCH "0")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

IF("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  MESSAGE(SEND_ERROR "In-source builds are not allowed, please build this in the bin directory!")
ENDIF("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")


# make vcflib
add_custom_command(
  OUTPUT ${PROJECT_BINARY_DIR}/vcflib.lib
  COMMAND make
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/externals/vcflib
  )

add_custom_target(
  vcflib_lib
  DEPENDS ${PROJECT_BINARY_DIR}/vcflib.lib
  )

# create library paths
set(VCFLIB_LIBRARIES
  ${PROJECT_SOURCE_DIR}/externals/vcflib/tabixpp/libtabix.a
  ${PROJECT_SOURCE_DIR}/externals/vcflib/libvcflib.a
  )

INCLUDE_DIRECTORIES(include src)


# add subfolders
ADD_SUBDIRECTORY(externals)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(tests)
