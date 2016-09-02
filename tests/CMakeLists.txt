# =================================
# graphite - test
#
# test/CMakeLists.txt
# =================================

# set header and source files
set(GRAPHITE_TEST_SOURCES
  gtest_main.cpp
)

# Where Google Test's .h files can be found.
include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/externals
  ${gtest_SOURCE_DIR}/include
  ${gtest_SOURCE_DIR}
  ${HTSLIB_INCLUDE_DIR}
  ${SAMTOOLS_INCLUDES}
  ${FASTAHACK_INCLUDE}
  ${ZLIB_INCLUDE_1}
  ${ZLIB_INCLUDE_2}
  ${CXXOPTS_INCLUDE}
 )



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()
set(CMAKE_BUILD_TYPE Debug)


#create binary
add_executable(graphite_tests
  ${GRAPHITE_TEST_SOURCES}
)

target_link_libraries(graphite_tests
    ${ZLIB_LIB}
    ${FASTAHACK_LIB}
	${HTSLIB_LIBRARY}
	gtest
	${BAMTOOLS_LIB}
	${BAMTOOLS_UTIL_LIB}
	${ZLIB_LIB}
	${VCFLIB_LIBRARIES}
	${GSSW_LIB}
	ssw_lib
    ${SCI_ZLIB_LIBRARY}
	${CORE_LIB}
	${FASTAHACK_LIB}
	${PLUGIN_LIBS}
	${SAMTOOLS_LIBRARIES}
	curl
	crypto
	${FASTAHACK_LIB}
	${ZLIB_LIB}
)

add_dependencies(graphite_tests ${GRAPHITE_EXTERNAL_PROJECT})
