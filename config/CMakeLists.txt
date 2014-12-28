
SET(SOMETHING "HELLO")

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/TestConfig.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/TestConfig.h"
)