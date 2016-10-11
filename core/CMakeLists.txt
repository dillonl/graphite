INCLUDE_DIRECTORIES(
    ${BAMTOOLS_INCLUDE}
    ${FASTAHACK_INCLUDE}
	${SCI_BOOST_INCLUDE}
	${HTSLIB_INCLUDE_DIR}
	${SAMTOOLS_INCLUDES}
	${TABIX_INCLUDE}
	${ZLIB_INCLUDE_1}
	${ZLIB_INCLUDE_2}
	${CXXOPTS_INCLUDE}
)

set(CMAKE_BUILD_TYPE Debug) # uncomment this when you need a debug build

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")

set(GRAPHITE_UTIL_SOURCES
  util/Params.cpp
  util/Utility.cpp
  )

set(GRAPHITE_SEQUENCE_SOURCES
  sequence/SequenceManager.cpp
  )

set(GRAPHITE_FILE_SOURCES
  file/BGZFFileWriter.cpp
  file/ASCIIFileWriter.cpp
  file/ASCIIFileReader.cpp
  file/ASCIIGZFileReader.cpp
  )

set(GRAPHITE_CORE_REFERENCE_SOURCES
  reference/Reference.cpp
  reference/FastaReference.cpp
  reference/FastaWriter.cpp
  )

set(GRAPHITE_CORE_REGION_SOURCES
  region/Region.cpp
  )

set(GRAPHITE_CORE_PATH_SOURCES
  path/Path.cpp
  path/PathManager.cpp
  )

set(GRAPHITE_CORE_MAPPING_SOURCES
  mapping/MappingManager.cpp
  )

set(GRAPHITE_CORE_ALIGNMENT_SOURCES
  alignment/AlignmentList.cpp
  alignment/BamAlignmentManager.cpp
  alignment/BamAlignmentReader.cpp
  alignment/SamtoolsAlignmentReader.cpp
  alignment/HTSLibAlignment.cpp
  alignment/HTSLibAlignmentReader.cpp
  alignment/HTSLibAlignmentManager.cpp
  )

set(GRAPHITE_CORE_VARIANT_SOURCES
  variant/Variant.cpp
  variant/VariantList.cpp
  variant/VCFManager.cpp
  variant/VCFFileReader.cpp
  variant/VCFHeader.cpp
  )

add_library(graphite_core STATIC
  ${GRAPHITE_UTIL_SOURCES}
  ${GRAPHITE_SEQUENCE_SOURCES}
  ${GRAPHITE_FILE_SOURCES}
  ${GRAPHITE_CORE_REFERENCE_SOURCES}
  ${GRAPHITE_CORE_REGION_SOURCES}
  ${GRAPHITE_CORE_ALIGNMENT_SOURCES}
  ${GRAPHITE_CORE_VARIANT_SOURCES}
  ${GRAPHITE_CORE_PATH_SOURCES}
  ${GRAPHITE_CORE_MAPPING_SOURCES}
  )

SET(CORE_LIB graphite_core CACHE STRING "THE CORE LIBRARY") #add graphite_core and graphite_utils to the list of core libs

add_dependencies(graphite_core ${TABIX_PROJECT} ${FASTAHACK_PROJECT} ${HTSLIB_PROJECT})
