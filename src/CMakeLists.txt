# =================================
# gwiz
#
# src/CMakeLists.txt
# =================================

# set header and source files
set(GWIZ_TOOLS_SOURCES
  tools/gwiz.cpp
)

set(GWIZ_UTILS_SOURCES
  core/utils/Parameters.cpp
  core/utils/file/ASCIIFileReader.cpp
)

set(GWIZ_CORE_SOURCES
  core/reference/Reference.cpp
  core/reference/FastaReference.cpp
  core/variants/Variant.cpp
  core/variants/VCFFileReader.cpp
  core/region/Region.cpp
  core/alignments/BamAlignment.cpp
  core/alignments/BamAlignmentReader.cpp
  core/alignments/BamAlignmentReaderManager.cpp
  core/alignments/BamAlignmentReaderPreload.cpp
  core/alignments/BamAlignmentReaderPreloadManager.cpp
)

set(GWIZ_GRAPH_SOURCES
)

include_directories(${CMAKE_SOURCE_DIR}/src/core/)
include_directories(${CMAKE_SOURCE_DIR}/externals/bamtools/src/)
include_directories(${CMAKE_SOURCE_DIR}/externals/)
include_directories(${SCI_BOOST_INCLUDE})

ADD_DEFINITIONS(-DBOOST_FALLTHROUGH)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

#set(CMAKE_BUILD_TYPE Debug)

add_library(gwiz_core ${GWIZ_CORE_SOURCES})
#add_library(gwiz_graph ${GWIZ_GRAPH_SOURCES})
add_library(gwiz_utils ${GWIZ_UTILS_SOURCES})

#create binary
add_executable(gwiz
  ${GWIZ_TOOLS_SOURCES}
)

# link libraries to executables
target_link_libraries(gwiz
	BamTools
	BamTools-utils
	${SCI_BOOST_LIBRARY}
	fastahack_lib
	ssw_lib
)

