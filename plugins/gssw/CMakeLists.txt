# =================================
# gssw
# =================================

set(GWIZ_GSSW_GRAPH_SOURCES
	graph/GSSWGraph.cpp
	graph/GraphManager.cpp
	graph/AlignmentReporter.cpp
	graph/AlignmentReport.cpp
	graph/GSSWAdjudicator.cpp
)

include_directories(${CMAKE_SOURCE_DIR}/src/)
include_directories(${CMAKE_SOURCE_DIR}/externals/bamtools/src/)
include_directories(${CMAKE_SOURCE_DIR}/externals/)
include_directories(${SCI_BOOST_INCLUDE})

#ADD_DEFINITIONS(-DBOOST_FALLTHROUGH)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")
set(CMAKE_BUILD_TYPE Debug)

add_library(gwiz_gssw_lib ${GWIZ_GSSW_GRAPH_SOURCES})

SET(TMP_PLUGINS ${PLUGIN_LIBS}) # store a temp plugin var
UNSET(PLUGIN_LIBS CACHE)  #unset the plugin_libs var
SET(PLUGIN_LIBS "${TMP_PLUGINS};gwiz_gssw_lib" CACHE LIST "A LIST OF THE PLUGIN LIBRARIES") #add gwiz_gssw_lib to the list of plugins

