# =================================
# adjudicator
# =================================

set(GRAPHITE_ADJUDICATOR_GRAPH_SOURCES
	graph/GSSWGraph.cpp
	graph/GraphManager.cpp
	graph/AlignmentReporter.cpp
	graph/AlignmentReport.cpp
	graph/GSSWAdjudicator.cpp
	graph/GSSWMapping.cpp
)

include_directories(${CMAKE_SOURCE_DIR}/src/)
include_directories(${CMAKE_SOURCE_DIR}/externals/bamtools/src/)
include_directories(${CMAKE_SOURCE_DIR}/externals/)
include_directories(${SCI_BOOST_INCLUDE} ${HTSLIB_INCLUDE_DIR} ${SAMTOOLS_INCLUDES})

#ADD_DEFINITIONS(-DBOOST_FALLTHROUGH)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")
#set(CMAKE_BUILD_TYPE Debug)

add_library(graphite_adjudicator_lib ${GRAPHITE_ADJUDICATOR_GRAPH_SOURCES})

SET(TMP_PLUGINS ${PLUGIN_LIBS}) # store a temp plugin var
UNSET(PLUGIN_LIBS CACHE)  #unset the plugin_libs var
SET(PLUGIN_LIBS "${TMP_PLUGINS};graphite_adjudicator_lib" CACHE LIST "A LIST OF THE PLUGIN LIBRARIES") #add graphite_adjudicator_lib to the list of plugins

