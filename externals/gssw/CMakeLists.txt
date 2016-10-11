# set header and source files
set(GSSW_SOURCES
  gssw.c
)

#set(CMAKE_BUILD_TYPE Debug) # uncomment this when you need a debug build

add_library(gssw_lib ${GSSW_SOURCES})

SET(GSSW_LIB gssw_lib CACHE STRING "THE GSSW LIBRARY")
