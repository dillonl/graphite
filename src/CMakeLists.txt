# =================================
# gwiz
#
# src/CMakeLists.txt
# =================================

set(GWIZ_UTILS_SOURCES
  utils/Parameters.cpp
  utils/file/ASCIIFileReader.cpp
  utils/file/ASCIIFileWriter.cpp
)

set(GWIZ_CORE_SOURCES
  reference/Reference.cpp
  reference/FastaReference.cpp
  variants/Variant.cpp
  variants/VCFFileReader.cpp
  variants/VariantListVCFPreloaded.cpp
  region/Region.cpp
  alignments/BamAlignment.cpp
  alignments/BamAlignmentReader.cpp
  alignments/BamAlignmentReaderManager.cpp
  alignments/BamAlignmentReaderPreload.cpp
  alignments/BamAlignmentReaderPreloadManager.cpp
)

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
add_library(gwiz_utils ${GWIZ_UTILS_SOURCES})
