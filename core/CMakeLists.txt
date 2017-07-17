INCLUDE_DIRECTORIES(
    ${BAMTOOLS_INCLUDE}
    ${FASTAHACK_INCLUDE}
	${CXXOPTS_INCLUDE}
	${TABIX_INCLUDE}
	${GSSW_INCLUDE}
	${ZLIB_INCLUDE}
	${CMAKE_CURRENT_SOURCE_DIR}/util
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
  util/gzstream.cpp
  )

set(GRAPHITE_SEQUENCE_SOURCES
  sequence/SequenceManager.cpp
  )

set(GRAPHITE_FILE_SOURCES
  file/BGZFFileWriter.cpp
  file/ASCIIFileWriter.cpp
  file/ASCIIFileReader.cpp
  file/ASCIIGZFileReader.cpp
  file/BamHeaderReader.cpp
  )

set(GRAPHITE_CORE_REFERENCE_SOURCES
  reference/Reference.cpp
  reference/FastaReference.cpp
  )

set(GRAPHITE_CORE_REGION_SOURCES
  region/Region.cpp
  )

set(GRAPHITE_CORE_MAPPING_SOURCES
  mapping/MappingManager.cpp
  mapping/GSSWMapping.cpp
  )

set(GRAPHITE_CORE_ADJUDICATOR_SOURCES
  adjudicator/GSSWAdjudicator.cpp
  )

set(GRAPHITE_CORE_ALIGNMENT_SOURCES
  alignment/AlignmentList.cpp
  alignment/BamAlignmentManager.cpp
  alignment/BamAlignmentReader.cpp
  alignment/AlignmentReporter.cpp
  alignment/AlignmentReport.cpp
  )

set(GRAPHITE_CORE_SAMPLE_SOURCES
  sample/Sample.cpp
  sample/SampleManager.cpp
  )

set(GRAPHITE_CORE_VARIANT_SOURCES
  variant/Variant.cpp
  variant/VariantList.cpp
  variant/VCFManager.cpp
  variant/VCFFileReader.cpp
  variant/VCFHeader.cpp
  )

set(GRAPHITE_CORE_GRAPH_SOURCES
  graph/ReferenceGraph.cpp
  graph/GSSWGraph.cpp
  graph/GraphManager.cpp
  )

add_library(graphite_core STATIC
  ${GRAPHITE_UTIL_SOURCES}
  ${GRAPHITE_SEQUENCE_SOURCES}
  ${GRAPHITE_FILE_SOURCES}
  ${GRAPHITE_CORE_REFERENCE_SOURCES}
  ${GRAPHITE_CORE_REGION_SOURCES}
  ${GRAPHITE_CORE_ALIGNMENT_SOURCES}
  ${GRAPHITE_CORE_SAMPLE_SOURCES}
  ${GRAPHITE_CORE_VARIANT_SOURCES}
  ${GRAPHITE_CORE_PATH_SOURCES}
  ${GRAPHITE_CORE_MAPPING_SOURCES}
  ${GRAPHITE_CORE_ADJUDICATOR_SOURCES}
  ${GRAPHITE_CORE_GRAPH_SOURCES}
  )

SET(CORE_LIB graphite_core CACHE STRING "THE CORE LIBRARY") #add graphite_core and graphite_utils to the list of core libs

TARGET_LINK_LIBRARIES(${CORE_LIB}
  ${BAMTOOLS_LIB}
  ${BAMTOOLS_UTIL_LIB}
  ${SCI_BOOST_LIBRARY}
  ${ZLIB_LIBRARY}
  ${FASTAHACK_LIB}
  ${TABIX_LIB}
  ${GSSW_LIB}
)

add_dependencies(${CORE_LIB} ${GRAPHITE_EXTERNAL_PROJECT})

