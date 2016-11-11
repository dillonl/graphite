# =================================
# graphite
#
# tools/CMakeLists.txt
# =================================

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)

# set header and source files
set(GRAPHITE_TOOLS_SOURCES
	graphite.cpp
)

# set header and source files
#set(PATHTRACE_TOOLS_SOURCES
#	pathtrace.cpp
#)

INCLUDE_DIRECTORIES(
  ${HTSLIB_INCLUDE_DIR}
  ${TABIX_INCLUDE}
  ${FASTAHACK_INCLUDE}
  ${ZLIB_INCLUDE_1}
  ${ZLIB_INCLUDE_2}
  ${CXXOPTS_INCLUDE}
  ${BAMTOOLS_INCLUDE}
)


#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -lstdc++fs -std=c++11") # testing out newer compiler
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

#set(CMAKE_BUILD_TYPE Debug) # uncomment this when you need a debug build

#create binary
add_executable(graphite
  ${GRAPHITE_TOOLS_SOURCES}
)
#add_executable(pathtrace
#  ${PATHTRACE_TOOLS_SOURCES}
#)

# link libraries to executables
target_link_libraries(graphite
  ${CURL_LIBRARIES}
  ${ZLIB_LIB}
  ${BAMTOOLS_LIB}
  ${BAMTOOLS_UTIL_LIB}
  ${FASTAHACK_LIB}
  ${ZLIB_LIB}
  ${HTSLIB_LIB}
  ${ZLIB_LIB}
  ${TABIX_LIB}
  ${PLUGIN_LIBS}
  ${CORE_LIB}
  ${FASTAHACK_LIB}
  ${ZLIB_LIB}
  ${ZLIB_LIB}
  ${GSSW_LIB}
  ${SCI_ZLIB_LIBRARY}
  ${SAMTOOLS_LIBRARIES}
  ${FASTAHACK_LIB}
  ${BAMTOOLS_LIB}
  ${BAMTOOLS_UTIL_LIB}
  ${ZLIB_LIB}
  ${CURL_LIBRARIES}
)

#target_link_libraries(pathtrace
#  ${ZLIB_LIB}
#  ${BAMTOOLS_LIB}
#  ${BAMTOOLS_UTIL_LIB}
#  ${ZLIB_LIB}
#  ${FASTAHACK_LIB}
#  ${TABIX_LIB}
#  ${PLUGIN_LIBS}
#  ${CORE_LIB}
#  ${FASTAHACK_LIB}
#  ${GSSW_LIB}
#  ${SCI_ZLIB_LIBRARY}
#  curl
#  crypto
#  ${TABIX_LIB}
#  ${FASTAHACK_LIB}
#  ${BAMTOOLS_LIB}
#  ${BAMTOOLS_UTIL_LIB}
#  ${ZLIB_LIB}
#)

add_dependencies(graphite ${GRAPHITE_EXTERNAL_PROJECT})
#add_dependencies(pathtrace ${GRAPHITE_EXTERNAL_PROJECT})


