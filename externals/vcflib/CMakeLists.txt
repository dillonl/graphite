# =================================
# gwiz
#
# src/CMakeLists.txt
# =================================

# set header and source files
set (VCF_LIB_SOURCES
  src/Variant.cpp
  src/split.cpp
)

add_library(vcf_lib SHARED VCF_LIB_SOURCES)

