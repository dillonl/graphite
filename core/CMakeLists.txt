
INCLUDE_DIRECTORIES(${SCI_BOOST_INCLUDE})

set(CMAKE_BUILD_TYPE Debug) # uncomment this when you need a debug build

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

set(GWIZ_UTIL_SOURCES
  util/Params.cpp
  )

set(GWIZ_SEQUENCE_SOURCES
  sequence/SequenceManager.cpp
  )

set(GWIZ_FILE_SOURCES
  file/ASCIIFileReader.cpp
  file/ASCIIGZFileReader.cpp
  )

set(GWIZ_CORE_REFERENCE_SOURCES
  reference/Reference.cpp
  reference/FastaReference.cpp
  )

set(GWIZ_CORE_REGION_SOURCES
  region/Region.cpp
  )

set(GWIZ_CORE_ALIGNMENT_SOURCES
  alignment/BamAlignment.cpp
  alignment/BamAlignmentReader.cpp
  alignment/BamAlignmentReaderManager.cpp
  alignment/BamAlignmentReaderPreload.cpp
  alignment/BamAlignmentReaderPreloadManager.cpp
  )

set(GWIZ_CORE_VARIANT_SOURCES
  variant/Variant.cpp
  variant/VariantList.cpp
  variant/VCFManager.cpp
  variant/VCFFileReader.cpp
  )

add_library(gwiz_core
  ${GWIZ_UTIL_SOURCES}
  ${GWIZ_SEQUENCE_SOURCES}
  ${GWIZ_FILE_SOURCES}
  ${GWIZ_CORE_REFERENCE_SOURCES}
  ${GWIZ_CORE_REGION_SOURCES}
  ${GWIZ_CORE_ALIGNMENT_SOURCES}
  ${GWIZ_CORE_VARIANT_SOURCES}
  )

SET(CORE_LIB gwiz_core CACHE STRING "THE CORE LIBRARY") #add gwiz_core and gwiz_utils to the list of core libs
