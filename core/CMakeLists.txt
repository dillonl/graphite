
INCLUDE_DIRECTORIES(${SCI_BOOST_INCLUDE})

set(CMAKE_BUILD_TYPE Debug) # uncomment this when you need a debug build

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

set(GWIZ_UTILS_SOURCES
  utils/Parameters.cpp
  utils/file/ASCIIFileReader.cpp
  utils/file/ASCIIGZFileReader.cpp
  )

set(GWIZ_CORE_REFERENCE_SOURCES
  reference/Reference.cpp
  reference/FastaReference.cpp
  )

set(GWIZ_CORE_REGION_SOURCES
  region/Region.cpp
  )

set(GWIZ_CORE_ALIGNMENTS_SOURCES
  alignments/BamAlignment.cpp
  alignments/BamAlignmentReader.cpp
  alignments/BamAlignmentReaderManager.cpp
  alignments/BamAlignmentReaderPreload.cpp
  alignments/BamAlignmentReaderPreloadManager.cpp
  )

set(GWIZ_CORE_VARIANTS_SOURCES
  variants/Variant.cpp
  variants/VariantListVCFPreloaded.cpp
  variants/VCFFileReader.cpp
  )

add_library(gwiz_core ${GWIZ_UTILS_SOURCES} ${GWIZ_CORE_REFERENCE_SOURCES} ${GWIZ_CORE_REGION_SOURCES} ${GWIZ_CORE_ALIGNMENTS_SOURCES} ${GWIZ_CORE_VARIANTS_SOURCES})

SET(CORE_LIB gwiz_core CACHE STRING "THE CORE LIBRARY") #add gwiz_core and gwiz_utils to the list of core libs
