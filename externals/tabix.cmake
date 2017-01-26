#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2015 Scientific Computing and Imaging Institute,
#  University of Utah.
# 
#  
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software. 
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.

SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

# Setting up external library for TABIX
SET(TABIX_PROJECT tabix_project CACHE INTERNAL "tabix project name")
SET(TABIX_DIR ${CMAKE_BINARY_DIR}/externals/tabix CACHE INTERNAL "tabix project directory")
ExternalProject_Add(${TABIX_PROJECT}
	GIT_REPOSITORY https://github.com/dillonl/tabix.git
	GIT_TAG 656f0c8ae8d8e94416c0a310cf22087a3880e585
	DEPENDS ${ZLIB_PROJECT}
	INSTALL_COMMAND ""
	UPDATE_COMMAND ""
	PREFIX ${TABIX_DIR}
    CMAKE_CACHE_ARGS
		-DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}
		-DZLIB_LIBRARY:PATH=${ZLIB_LIBRARY}
		-DZLIB_INCLUDE:PATH=${ZLIB_INCLUDE}
)


ExternalProject_Get_Property(${TABIX_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${TABIX_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${TABIX_PROJECT} BINARY_DIR)

SET(TABIX_LIB ${BINARY_DIR}/libtabix.a CACHE INTERNAL "TABIX Library")
SET(TABIX_INCLUDE ${SOURCE_DIR} CACHE INTERNAL "TABIX Include")

