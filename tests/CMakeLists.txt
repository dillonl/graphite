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
  ${GSSW_INCLUDE}
  ${ZLIB_INCLUDE}
  ${BAMTOOLS_INCLUDE}
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/externals
  ${CMAKE_SOURCE_DIR}/core/util
  ${GTEST_SOURCE_DIR}/include
  ${GTEST_SOURCE_DIR}
  ${GMOCK_SOURCE_DIR}
  ${FASTAHACK_INCLUDE}
  ${SCI_BOOST_INCLUDE}
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
  ${CORE_LIB}
  ${GRAPHITE_ADJUDICATOR}
  ${GTEST_LIB}
  ${GMOCK_LIB}
)

add_dependencies(graphite_tests ${GRAPHITE_EXTERNAL_PROJECT})
