# =================================
# graphite
#
# tools/CMakeLists.txt
# =================================

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)

# set header and source files
set(GRAPHITE_TOOLS_SOURCES
	graphite.cpp
)

# set header and source files
set(PATHTRACE_TOOLS_SOURCES
	pathtrace.cpp
)

ADD_DEFINITIONS(-DBOOST_FALLTHROUGH)
INCLUDE_DIRECTORIES(${SCI_BOOST_INCLUDE})


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

set(CMAKE_BUILD_TYPE Debug) # uncomment this when you need a debug build

#create binary
add_executable(graphite
  ${GRAPHITE_TOOLS_SOURCES}
)
add_executable(pathtrace
  ${PATHTRACE_TOOLS_SOURCES}
)

# link libraries to executables
target_link_libraries(graphite
  BamTools
  BamTools-utils
  ${PLUGIN_LIBS}
  ${CORE_LIB}
  fastahack_lib
  ${SCI_BOOST_LIBRARY}
  ${GSSW_LIB}
)

target_link_libraries(pathtrace
  BamTools
  BamTools-utils
  ${PLUGIN_LIBS}
  ${CORE_LIB}
  fastahack_lib
  ${SCI_BOOST_LIBRARY}
  ${GSSW_LIB}
)

