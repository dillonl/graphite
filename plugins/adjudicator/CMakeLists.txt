# =================================
# vg
# =================================

set(GRAPHITE_ADJUDICATOR_GRAPH_SOURCES
	graph/AdjudicatorGraph.cpp
	graph/VariantContig.cpp
)

include_directories(${CMAKE_SOURCE_DIR}/src/)
include_directories(${CMAKE_SOURCE_DIR}/externals/bamtools/src/)
include_directories(${CMAKE_SOURCE_DIR}/externals/)
include_directories(${SCI_BOOST_INCLUDE})

#ADD_DEFINITIONS(-DBOOST_FALLTHROUGH)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

#set(CMAKE_BUILD_TYPE Debug)

add_library(adjudicator_lib ${GRAPHITE_ADJUDICATOR_GRAPH_SOURCES})

SET(TMP_PLUGINS ${PLUGIN_LIBS}) # store a temp plugin var
UNSET(PLUGIN_LIBS CACHE)  #unset the plugin_libs var
SET(PLUGIN_LIBS "${TMP_PLUGINS};adjudicator_lib" CACHE LIST "A LIST OF THE PLUGIN LIBRARIES") #add adjudicator_lib to the list of plugins

#create binary
#add_executable(graphite
#  ${GRAPHITE_TOOLS_SOURCES}
#)

# link libraries to executables
#target_link_libraries(graphite BamTools BamTools-utils ${SCI_BOOST_LIBRARY})

