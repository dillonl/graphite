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

INCLUDE_DIRECTORIES(include config externals plugins ${PROJECT_SOURCE_DIR})

UNSET(CORE_LIBS CACHE) #unset this each time
SET(CORE_LIBS CACHE LIST "A LIST OF THE PLUGIN LIBRARIES")

# add subfolders
ADD_SUBDIRECTORY(externals)
ADD_SUBDIRECTORY(config)
ADD_SUBDIRECTORY(core)
ADD_SUBDIRECTORY(plugins)
ADD_SUBDIRECTORY(tools)
ADD_SUBDIRECTORY(tests)
