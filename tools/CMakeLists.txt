# =================================
# gwiz
#
# tools/CMakeLists.txt
# =================================

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)

# set header and source files
set(GWIZ_TOOLS_SOURCES
  gssw.cpp
)

ADD_DEFINITIONS(-DBOOST_FALLTHROUGH)
INCLUDE_DIRECTORIES(${SCI_BOOST_INCLUDE})


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

set(CMAKE_BUILD_TYPE Debug) # uncomment this when you need a debug build

#create binary
add_executable(gssw
  ${GWIZ_TOOLS_SOURCES}
)

# link libraries to executables
target_link_libraries(gssw
  BamTools
  BamTools-utils
  ${CORE_LIB}
  fastahack_lib
  ${SCI_BOOST_LIBRARY}
  ${PLUGIN_LIBS}
  ${GSSW_LIB}
)

