# =================================
# graphite - test
#
# test/CMakeLists.txt
# =================================

# set header and source files
set(GRAPHITE_TEST_SOURCES
  gtest_main.cpp
)

message("samtools: ${SAMTOOLS_INCLUDES}")

# Where Google Test's .h files can be found.
include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/externals
  ${gtest_SOURCE_DIR}/include
  ${gtest_SOURCE_DIR}
  ${SCI_BOOST_INCLUDE}
  ${HTSLIB_INCLUDE_DIR}
  ${SAMTOOLS_INCLUDES}
 )



set(Boost_USE_MULTITHREADED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()
#set(CMAKE_BUILD_TYPE Debug)


#create binary
add_executable(graphite_tests
  ${GRAPHITE_TEST_SOURCES}
)

target_link_libraries(graphite_tests
	${HTSLIB_LIBRARY}
	gtest
	BamTools
	BamTools-utils
	${VCFLIB_LIBRARIES}
	fastahack_lib
	${GSSW_LIB}
	ssw_lib
	${SCI_BOOST_LIBRARY}
    ${SCI_ZLIB_LIBRARY}
	${CORE_LIB}
	${PLUGIN_LIBS}
	${SAMTOOLS_LIBRARIES}
)

add_dependencies(graphite_tests ${SAMTOOLS_PROJECT} ${HTSLIB_PROJECT})
