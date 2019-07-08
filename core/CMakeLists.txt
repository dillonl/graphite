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
  util/GraphPrinter.cpp
  util/Params.cpp
  util/Utility.cpp
  util/gzstream.cpp
  )

set(GRAPHITE_CORE_REFERENCE_SOURCES
  reference/FastaReference.cpp
  )

set(GRAPHITE_CORE_REGION_SOURCES
  region/Region.cpp
  )

set(GRAPHITE_CORE_VCF_SOURCES
  vcf/VCFReader.cpp
  vcf/VCFWriter.cpp
  vcf/Variant.cpp
  )

set(GRAPHITE_CORE_BAM_SOURCES
  bam/BamReader.cpp
  )

set(GRAPHITE_CORE_GRAPH_PROCESSOR_SOURCES
  graph/GraphProcessor.cpp
  graph/ReferenceGraph.cpp
  graph/Graph.cpp
  graph/Traceback.cpp
  graph/Node.cpp
  )

set(GRAPHITE_CORE_SAMPLE_SOURCES
  sample/Sample.cpp
  )

set(GRAPHITE_CORE_ALLELE_SOURCES
  allele/Allele.cpp
  )

add_library(graphite_core STATIC
  ${GRAPHITE_UTIL_SOURCES}
  ${GRAPHITE_CORE_REFERENCE_SOURCES}
  ${GRAPHITE_CORE_REGION_SOURCES}
  ${GRAPHITE_CORE_VCF_SOURCES}
  ${GRAPHITE_CORE_BAM_SOURCES}
  ${GRAPHITE_CORE_GRAPH_PROCESSOR_SOURCES}
  ${GRAPHITE_CORE_SAMPLE_SOURCES}
  ${GRAPHITE_CORE_ALLELE_SOURCES}
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

