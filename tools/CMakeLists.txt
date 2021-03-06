# =================================
# graphite
#
# tools/CMakeLists.txt
# =================================

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)


set(GRAPHITE_TOOLS_SOURCES
	graphite.cpp
)

# set header and source files
#set(PATHTRACE_TOOLS_SOURCES
#	pathtrace.cpp
#)

ADD_DEFINITIONS(-DBOOST_FALLTHROUGH)
INCLUDE_DIRECTORIES(
  ${ZLIB_INCLUDE}
  ${TABIX_INCLUDE}
  ${GSSW_INCLUDE}
  ${FASTAHACK_INCLUDE}
  ${CXXOPTS_INCLUDE}
  ${HTSLIB_INCLUDE}
)


#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -lstdc++fs -std=c++11") # testing out newer compiler
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

set(CMAKE_BUILD_TYPE Debug) # uncomment this when you need a debug build

add_executable(graphite
  ${GRAPHITE_TOOLS_SOURCES}
)

# link libraries to executables
target_link_libraries(graphite
  ${CORE_LIB}
  ${GRAPHITE_ADJUDICATOR}
)

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
target_link_libraries(graphite ${PYTHON_LIBRARIES})

add_dependencies(graphite ${GRAPHITE_EXTERNAL_PROJECT})
#install(TARGETS graphite DESTINATION bin)