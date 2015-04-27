#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2009 Scientific Computing and Imaging Institute,
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
#


# THIS IS A SELECTION OF THE BOOST LIBRARY
# IT HAS BEEN STRIPPED DOWN TO ITS BARE ESSENTIALS
# THIS CMAKEFILE GENERATES A SIMPLE BOOST LIBRARY THAT
# IS STATICALLY LINKED.



###########################################
# Boost project

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(Boost)

###########################################
# Source files for zlib

SET(BOOST_SRCS
  libs/chrono/src/chrono.cpp
  libs/chrono/src/process_cpu_clocks.cpp
  libs/chrono/src/thread_clock.cpp

  libs/date_time/src/gregorian/date_generators.cpp 
  libs/date_time/src/gregorian/greg_month.cpp 
  libs/date_time/src/gregorian/greg_weekday.cpp 
  libs/date_time/src/gregorian/gregorian_types.cpp 
  libs/date_time/src/posix_time/posix_time_types.cpp 

  libs/filesystem/src/codecvt_error_category.cpp
  libs/filesystem/src/operations.cpp
  libs/filesystem/src/path.cpp
  libs/filesystem/src/path_traits.cpp
  libs/filesystem/src/portability.cpp
  libs/filesystem/src/unique_path.cpp
  libs/filesystem/src/utf8_codecvt_facet.cpp
  libs/filesystem/src/windows_file_codecvt.cpp
  libs/filesystem/src/windows_file_codecvt.hpp

#  libs/detail/utf8_codecvt_facet.cpp
  
  libs/graph/src/graphml.cpp
  libs/graph/src/read_graphviz_new.cpp
#  libs/graph_parallel/src/mpi_process_group.cpp
#  libs/graph_parallel/src/tag_allocator.cpp
  
  libs/iostreams/src/file_descriptor.cpp 
  libs/iostreams/src/mapped_file.cpp
# the PREFIX mangling does do not combine with boost
# need to alter code to fix this
  libs/iostreams/src/zlib.cpp
  libs/iostreams/src/bzip2.cpp
  libs/iostreams/src/gzip.cpp

  libs/program_options/src/cmdline.cpp 
  libs/program_options/src/config_file.cpp 
  libs/program_options/src/convert.cpp 
  libs/program_options/src/options_description.cpp 
  libs/program_options/src/parsers.cpp 
  libs/program_options/src/positional_options.cpp 
  libs/program_options/src/split.cpp 
  libs/program_options/src/utf8_codecvt_facet.cpp 
  libs/program_options/src/value_semantic.cpp 
  libs/program_options/src/variables_map.cpp 
  libs/program_options/src/winmain.cpp
  
  libs/regex/src/c_regex_traits.cpp
  libs/regex/src/cpp_regex_traits.cpp
  libs/regex/src/cregex.cpp
  libs/regex/src/fileiter.cpp
  libs/regex/src/icu.cpp
  libs/regex/src/instances.cpp
  libs/regex/src/posix_api.cpp
  libs/regex/src/regex.cpp
  libs/regex/src/regex_debug.cpp
  libs/regex/src/regex_raw_buffer.cpp
  libs/regex/src/regex_traits_defaults.cpp
  libs/regex/src/static_mutex.cpp
  libs/regex/src/w32_regex_traits.cpp
  libs/regex/src/wc_regex_traits.cpp
  libs/regex/src/wide_posix_api.cpp
  libs/regex/src/winstances.cpp 
  libs/regex/src/usinstances.cpp   

#  libs/serialization/src/archive_exception.cpp  
#  libs/serialization/src/basic_archive.cpp
#  libs/serialization/src/basic_iarchive.cpp
#  libs/serialization/src/basic_iserializer.cpp
#  libs/serialization/src/basic_oarchive.cpp
#  libs/serialization/src/basic_oserializer.cpp
#  libs/serialization/src/basic_pointer_iserializer.cpp
#  libs/serialization/src/basic_pointer_oserializer.cpp
#  libs/serialization/src/basic_serializer_map.cpp
#  libs/serialization/src/basic_text_iprimitive.cpp
#  libs/serialization/src/basic_text_oprimitive.cpp
#  libs/serialization/src/basic_text_wiprimitive.cpp
#  libs/serialization/src/basic_text_woprimitive.cpp
#  libs/serialization/src/basic_xml_archive.cpp
#  libs/serialization/src/binary_iarchive.cpp
#  libs/serialization/src/binary_oarchive.cpp
#  libs/serialization/src/binary_wiarchive.cpp
#  libs/serialization/src/binary_woarchive.cpp
#  libs/serialization/src/codecvt_null.cpp
#  libs/serialization/src/extended_type_info.cpp
#  libs/serialization/src/extended_type_info_typeid.cpp
#  libs/serialization/src/extended_type_info_no_rtti.cpp
#  libs/serialization/src/polymorphic_iarchive.cpp
#  libs/serialization/src/polymorphic_oarchive.cpp
#  libs/serialization/src/shared_ptr_helper.cpp
#  libs/serialization/src/stl_port.cpp
#  libs/serialization/src/text_iarchive.cpp
#  libs/serialization/src/text_oarchive.cpp
#  libs/serialization/src/text_wiarchive.cpp
#  libs/serialization/src/text_woarchive.cpp
#  libs/serialization/src/utf8_codecvt_facet.cpp
#  libs/serialization/src/void_cast.cpp
#  libs/serialization/src/xml_archive_exception.cpp
#  libs/serialization/src/xml_grammar.cpp
#  libs/serialization/src/xml_iarchive.cpp
#  libs/serialization/src/xml_oarchive.cpp
#  libs/serialization/src/xml_wgrammar.cpp
#  libs/serialization/src/xml_wiarchive.cpp
#  libs/serialization/src/xml_woarchive.cpp

  libs/system/src/error_code.cpp  
)

# Process library not needed for SCIRun and Seg3D, so force off for now
SET(BUILD_WITH_PROCESS OFF CACHE INTERNAL "Optionally build Boost process (unapproved project from sandbox) library." FORCE)
IF(BUILD_WITH_PROCESS)
  SET(BOOST_SRCS ${BOOST_SRCS}
    libs/processes/src/child.cpp
    libs/processes/src/environment.cpp
    libs/processes/src/launch_shell.cpp
    libs/processes/src/launch.cpp
    libs/processes/src/pipe_end.cpp
    libs/processes/src/pipe.cpp
    libs/processes/src/posix_util.cpp
    libs/processes/src/utility.cpp
    libs/processes/src/windows_util.cpp
  )
ENDIF(BUILD_WITH_PROCESS)

IF(BUILD_WITH_PYTHON)
  SET(BOOST_SRCS ${BOOST_SRCS}
   libs/python/src/numeric.cpp
   libs/python/src/list.cpp
   libs/python/src/long.cpp
   libs/python/src/dict.cpp
   libs/python/src/tuple.cpp
   libs/python/src/str.cpp
   libs/python/src/slice.cpp

   libs/python/src/converter/from_python.cpp
   libs/python/src/converter/registry.cpp
   libs/python/src/converter/type_id.cpp
   libs/python/src/object/enum.cpp
   libs/python/src/object/class.cpp
   libs/python/src/object/function.cpp
   libs/python/src/object/inheritance.cpp
   libs/python/src/object/life_support.cpp
   libs/python/src/object/pickle_support.cpp
   libs/python/src/errors.cpp
   libs/python/src/module.cpp
   libs/python/src/converter/builtin_converters.cpp
   libs/python/src/converter/arg_to_python_base.cpp
   libs/python/src/object/iterator.cpp
   libs/python/src/object/stl_iterator.cpp
   libs/python/src/object_protocol.cpp
   libs/python/src/object_operators.cpp
   libs/python/src/wrapper.cpp
   libs/python/src/import.cpp
   libs/python/src/exec.cpp
   libs/python/src/object/function_doc_signature.cpp
  )
  ADD_DEFINITIONS(-DBOOST_PYTHON_STATIC_LIB -DBOOST_PYTHON_SOURCE)
  INCLUDE_DIRECTORIES(${SCI_PYTHON_INCLUDE})
ENDIF(BUILD_WITH_PYTHON)

IF (WIN32)
  SET(BOOST_SRCS
   ${BOOST_SRCS}
   libs/thread/src/win32/thread.cpp
   libs/thread/src/win32/tss_dll.cpp
   libs/thread/src/win32/tss_pe.cpp)
ELSE (WIN32)
  SET(BOOST_SRCS
   ${BOOST_SRCS}
   libs/thread/src/pthread/thread.cpp 
   libs/thread/src/pthread/once.cpp)
ENDIF (WIN32)


###########################################
# Set the name of the library

SET(SCI_BOOST_LIBRARY boost CACHE INTERNAL "Boost library name." FORCE)

###########################################
# Include directory of zlib

SET(SCI_BOOST_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR} 
    CACHE INTERNAL "Boost include directories." FORCE)
INCLUDE_DIRECTORIES(${SCI_BOOST_INCLUDE})
INCLUDE_DIRECTORIES(${SCI_CRYPTO_INCLUDE})
INCLUDE_DIRECTORIES(${SCI_CRYPTO_INCLUDE}/../include)

###########################################
# Defines copied from the CMakefiles provided
# by boost

ADD_DEFINITIONS(-DDATE_TIME_INLINE)
ADD_DEFINITIONS(-DBOOST_DATE_TIME_STATIC_LINK)
ADD_DEFINITIONS(-DBOOST_FILESYSTEM_STATIC_LINK)
ADD_DEFINITIONS(-DBOOST_THREAD_BUILD_LIB)
ADD_DEFINITIONS(-DZ_PREFIX)
ADD_DEFINITIONS(-DBOOST_ALL_NO_LIB)

###########################################
# Add the library

INCLUDE_DIRECTORIES(${SCI_ZLIB_INCLUDE})
ADD_LIBRARY(${SCI_BOOST_LIBRARY} STATIC ${BOOST_SRCS})
TARGET_LINK_LIBRARIES(${SCI_BOOST_LIBRARY} ${SCI_ZLIB_LIBRARY})
IF(BUILD_WITH_PYTHON)
  TARGET_LINK_LIBRARIES(${SCI_BOOST_LIBRARY} ${SCI_PYTHON_LIBRARY})
  SET_TARGET_PROPERTIES(${SCI_BOOST_LIBRARY}
                        PROPERTIES
                        COMPILE_DEFINITIONS_DEBUG "BOOST_DEBUG_PYTHON")
ENDIF(BUILD_WITH_PYTHON)

###########################################
# Add iostreams to the boost library
#TARGET_LINK_LIBRARIES(${SCI_BOOST_LIBRARY} ${SCI_ZLIB_LIBRARY})



# required by threading library for clock_gettime
IF(UNIX AND NOT APPLE AND NOT CYGWIN)
  TARGET_LINK_LIBRARIES(${SCI_BOOST_LIBRARY} rt)
ENDIF(UNIX AND NOT APPLE AND NOT CYGWIN)


###########################################
# Enable PIC if necessary

IF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)
  SET(ENABLE_PIC ON)
ENDIF(BUILD_SHARED_LIBS AND UNIX AND NOT APPLE AND NOT CYGWIN)

IF(ENABLE_PIC)
  SET_TARGET_PROPERTIES(${SCI_BOOST_LIBRARY}
                        PROPERTIES
                          COMPILE_FLAGS "-fPIC"
                       )
ENDIF(ENABLE_PIC)

###########################################
# Set Windows options (Visual Studio)
# Force Visual Studio to link against the 
# version that is in the redistributable
# Options for Visual Studio to limit 
# the number of warnings.

IF(WIN32)
  ADD_DEFINITIONS(-D_BIND_TO_CURRENT_VCLIBS_VERSION=1)

  IF(NOT UNIX AND NOT BORLAND AND NOT MINGW)
    SET_TARGET_PROPERTIES(${SCI_BOOST_LIBRARY}
                            PROPERTIES
                            COMPILE_FLAGS "-W3 -wd4506 -wd4800 -wd4244 -wd4267 -wd4005"
                            LINK_FLAGS "/DYNAMICBASE:NO")
  ENDIF(NOT UNIX AND NOT BORLAND AND NOT MINGW)
ENDIF(WIN32)
