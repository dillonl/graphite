# =================================
# vg
# =================================

set(GWIZ_VG_GRAPH_SOURCES
	graph/VariantGraph.cpp
)

include_directories(${CMAKE_SOURCE_DIR}/src/)
include_directories(${CMAKE_SOURCE_DIR}/externals/bamtools/src/)
include_directories(${CMAKE_SOURCE_DIR}/externals/)
include_directories(${SCI_BOOST_INCLUDE})

#ADD_DEFINITIONS(-DBOOST_FALLTHROUGH)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")
set(CMAKE_BUILD_TYPE Debug)

add_library(vg_lib ${GWIZ_VG_GRAPH_SOURCES})

SET(TMP_PLUGINS ${PLUGIN_LIBS}) # store a temp plugin var
UNSET(PLUGIN_LIBS CACHE)  #unset the plugin_libs var
SET(PLUGIN_LIBS "${TMP_PLUGINS};vg_lib" CACHE LIST "A LIST OF THE PLUGIN LIBRARIES") #add vg_lib to the list of plugins

#create binary
#add_executable(gwiz
#  ${GWIZ_TOOLS_SOURCES}
#)

# link libraries to executables
#target_link_libraries(gwiz BamTools BamTools-utils ${SCI_BOOST_LIBRARY})

