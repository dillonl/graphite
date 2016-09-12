# build htslib
SET(HTSLIB_DIR ${CMAKE_BINARY_DIR}/externals/htslib CACHE INTERNAL "htslib project directory")
SET(HTSLIB_PROJECT htslib_project CACHE INTERNAL "htslib project name")
ExternalProject_Add(${HTSLIB_PROJECT}
	DEPENDS ${ZLIB_PROJECT}
    PREFIX ${HTSLIB_DIR}
    GIT_REPOSITORY "https://github.com/samtools/htslib.git"
    GIT_TAG 11661a57306a048528d0a223c86c62bcdc20eb18
    BUILD_IN_SOURCE 1
	UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND  "make"
    INSTALL_COMMAND  ""
    LOG_DOWNLOAD 0
	LOG_UPDATE 0
    LOG_CONFIGURE 0
    LOG_BUILD 0
    LOG_TEST 0
    LOG_INSTALL 0
	CMAKE_CACHE_ARGS
        -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}
)

#target_link_libraries (htslib curl crypto)

#include_directories(${HTSLIB_DIR}/src/HTSLIB_PROJECT/)

ExternalProject_Get_Property(${HTSLIB_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${HTSLIB_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${HTSLIB_PROJECT} BINARY_DIR)

SET(HTSLIB_INCLUDE_DIR ${SOURCE_DIR} CACHE INTERNAL "htslib include")
SET(HTSLIB_LIB ${SOURCE_DIR}/libhts.a CACHE INTERNAL "htslib Library")

