# =================================
# gwiz
#
# src/CMakeLists.txt
# =================================

message ("building vcflib")
execute_process(COMMAND make -j10
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# set header and source files
#set (VCF_LIB_SOURCES
#  src/Variant.cpp
#  src/split.cpp
#)

#add_library(vcf_lib SHARED VCF_LIB_SOURCES)

