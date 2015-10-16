# =================================
# pathtrace
# =================================

set(GRAPHITE_PATHTRACE_GRAPH_SOURCES
	graph/VariantGraph.cpp
)

include_directories(${CMAKE_SOURCE_DIR}/src/)
include_directories(${CMAKE_SOURCE_DIR}/externals/bamtools/src/)
include_directories(${CMAKE_SOURCE_DIR}/externals/)
include_directories(${SCI_BOOST_INCLUDE})

#ADD_DEFINITIONS(-DBOOST_FALLTHROUGH)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")
#set(CMAKE_BUILD_TYPE Debug)

add_library(pathtrace_lib ${GRAPHITE_PATHTRACE_GRAPH_SOURCES})

SET(TMP_PLUGINS ${PLUGIN_LIBS}) # store a temp plugin var
UNSET(PLUGIN_LIBS CACHE)  #unset the plugin_libs var
SET(PLUGIN_LIBS "${TMP_PLUGINS};pathtrace_lib" CACHE LIST "A LIST OF THE PLUGIN LIBRARIES") #add pathtrace_lib to the list of plugins

#create binary
#add_executable(graphite
#  ${GRAPHITE_TOOLS_SOURCES}
#)

# link libraries to executables
#target_link_libraries(graphite BamTools BamTools-utils ${SCI_BOOST_LIBRARY})

