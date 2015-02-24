# =================================
# gwiz - test
#
# text/CMakeLists.txt
# =================================

# set header and source files
set(GWIZ_TEST_SOURCES
  gtest_main.cpp
)

# Where Google Test's .h files can be found.
include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/externals
  ${gtest_SOURCE_DIR}/include
  ${gtest_SOURCE_DIR}
  ${SCI_BOOST_INCLUDE})

set(Boost_USE_MULTITHREADED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()
#set(CMAKE_BUILD_TYPE Debug)


#create binary
add_executable(gwiz_tests
  ${GWIZ_TEST_SOURCES}
)

MESSAGE("PLUGIN: ${PLUGIN_LIBS}")
target_link_libraries(gwiz_tests
	gtest
	BamTools
	BamTools-utils
	${VCFLIB_LIBRARIES}
	gwiz_core
	gwiz_utils
	${SCI_BOOST_LIBRARY}
	${PLUGIN_LIBS}
	fastahack_lib
	gssw_lib
	ssw_lib
)
