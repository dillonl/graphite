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

# Setting up external library gssw, we don't build it because we only need the include directories

SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

SET(GSSW_PROJECT gssw_project CACHE INTERNAL "gssw project name")
SET(GSSW_DIR ${CMAKE_BINARY_DIR}/externals/gssw CACHE INTERNAL "gssw project directory")
SET(GSSW_LIB)
ExternalProject_Add(${GSSW_PROJECT}
	GIT_REPOSITORY https://github.com/dillonl/gssw-1.git
	GIT_TAG db7e38552538ce703c98cee803e79df84cd2b148 #lock in the commit id so we don't this doesn't break in the future
	CONFIGURE_COMMAND ""
	BUILD_COMMAND "make"
	INSTALL_COMMAND ""
	UPDATE_COMMAND ""
	BUILD_IN_SOURCE 1
	PREFIX ${GSSW_DIR}
    LOG_DOWNLOAD 0
	LOG_UPDATE 0
    LOG_CONFIGURE 0
    LOG_BUILD 0
    LOG_TEST 0
    LOG_INSTALL 0
)

ExternalProject_Get_Property(${GSSW_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${GSSW_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${GSSW_PROJECT} BINARY_DIR)

SET(GSSW_LIB ${BINARY_DIR}/lib/libgssw.a CACHE INTERNAL "GSSW Lib")
SET(GSSW_INCLUDE ${SOURCE_DIR}/src/ CACHE INTERNAL "GSSW Include")

