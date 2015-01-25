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

INCLUDE_DIRECTORIES(include src config plugins)

SET(PLUGIN_LIBS)

# add subfolders
ADD_SUBDIRECTORY(config)
ADD_SUBDIRECTORY(externals)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(plugins)
ADD_SUBDIRECTORY(tests)
