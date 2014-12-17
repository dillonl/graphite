# =================================
# gwiz
#
# src/CMakeLists.txt
# =================================

# set header and source files
set(GWIZ_SOURCES
  utils/Parameters.cpp
  tools/gwiz.cpp
)

set(GWIZ_CORE_SOURCES
  core/Reference.cpp
  core/Variant.cpp
)

set(GWIZ_GRAPH_SOURCES
  graph/IGraph.cpp
  graph/SparseGraph.cpp
)

include_directories(${CMAKE_SOURCE_DIR}/externals/bamtools/src/)
include_directories(${CMAKE_SOURCE_DIR}/externals/vcflib)
include_directories(${CMAKE_SOURCE_DIR}/externals/vcflib/src)
include_directories(${CMAKE_SOURCE_DIR}/externals/)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")
set(CMAKE_BUILD_TYPE Debug)

add_library(gwiz_core ${GWIZ_CORE_SOURCES})
add_library(gwiz_graph ${GWIZ_GRAPH_SOURCES})

#create binary
add_executable(gwiz
  ${GWIZ_SOURCES}
)

# link libraries to executables
target_link_libraries(gwiz BamTools BamTools-utils ${VCFLIB_LIBRARIES})

