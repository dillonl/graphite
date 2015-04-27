#   zlib source code license
#
#   Copyright (C) 1995-2005 Jean-loup Gailly and Mark Adler
# 
#   This software is provided 'as-is', without any express or implied
#   warranty.  In no event will the authors be held liable for any damages
#   arising from the use of this software.
# 
#   Permission is granted to anyone to use this software for any purpose,
#   including commercial applications, and to alter it and redistribute it
#   freely, subject to the following restrictions:
# 
#   1. The origin of this software must not be misrepresented; you must not
#      claim that you wrote the original software. If you use this software
#      in a product, an acknowledgment in the product documentation would be
#      appreciated but is not required.
#   2. Altered source versions must be plainly marked as such, and must not be
#      misrepresented as being the original software.
#   3. This notice may not be removed or altered from any source distribution.
# 
#   Jean-loup Gailly        Mark Adler
#   jloup@gzip.org          madler@alumni.caltech.edu
# 
# 
#   The data format used by the zlib library is described by RFCs (Request for
#   Comments) 1950 to 1952 in the files http://www.ietf.org/rfc/rfc1950.txt
#   (zlib format), rfc1951.txt (deflate format) and rfc1952.txt (gzip format).

# CMake build scripts provided by Scientific Computing and Imaging Institute
# Copyright (c) 2009 Scientific Computing and Imaging Institute, University of Utah.

###########################################
# ZLib project

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(ZLIB)

###########################################
# Test for some required system information.

INCLUDE (CMakeBackwardCompatibilityC)

###########################################
# Source files for zlib

SET(ZLIB_SRCS
      adler32.c
      compress.c
      crc32.c
      deflate.c
      gzio.c
      infback.c
      inffast.c
      inflate.c
      inftrees.c
      trees.c
      uncompr.c
      zutil.c
    )

###########################################
# for windows add the .def and 
# .rc files to the source list
# if building shared libs

IF(WIN32)
  IF(BUILD_SHARED_LIBS)
    IF(NOT UNIX AND NOT BORLAND AND NOT MINGW)
      SET(ZLIB_SRCS ${ZLIB_SRCS} win32/zlib.def win32/zlib1.rc)
    ENDIF(NOT UNIX AND NOT BORLAND AND NOT MINGW)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(WIN32)

###########################################
# Set the name of the library

IF(WIN32)
  SET(SCI_ZLIB_LIBRARY zlib1 CACHE INTERNAL "Zlib library name." FORCE)
ELSE(WIN32)
  SET(SCI_ZLIB_LIBRARY zlib CACHE INTERNAL "Zlib library name." FORCE)
ENDIF(WIN32)

###########################################
# Include directory of zlib

SET(SCI_ZLIB_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR} 
    CACHE INTERNAL "Zlib include directories." FORCE)
INCLUDE_DIRECTORIES(${SCI_ZLIB_INCLUDE})

###########################################
# Force Zlib names to be mangled to enusre
# they do not conflict with other packages
# such as QT

ADD_DEFINITIONS(-DZ_PREFIX)

SET(SCI_ZLIB_MANGLE 1 CACHE INTERNAL "Tell the system that we have a mangled zlib." FORCE)

###########################################
# Add the library

ADD_LIBRARY(${SCI_ZLIB_LIBRARY} STATIC ${ZLIB_SRCS})


###########################################
# Set Windows options (Visual Studio)
# Force Visual Studio to link against the 
# version that is in the redistributable
# Options for Visual Studio to limit 
# the number of warnings.

IF(WIN32)
  ADD_DEFINITIONS(-D_BIND_TO_CURRENT_VCLIBS_VERSION=1)

  IF(NOT UNIX AND NOT BORLAND AND NOT MINGW)
    SET_TARGET_PROPERTIES(${SCI_ZLIB_LIBRARY}
                            PROPERTIES
                            COMPILE_FLAGS "/Ob1 /Gy /GF -W3 -wd4267 -wd4244"
                            LINK_FLAGS "/DYNAMICBASE:NO")
  ENDIF(NOT UNIX AND NOT BORLAND AND NOT MINGW)
ENDIF(WIN32)

IF(CMAKE_COMPILER_IS_GNUCC)
  SET_TARGET_PROPERTIES(${SCI_ZLIB_LIBRARY}
                        PROPERTIES
                        COMPILE_FLAGS "-w")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
