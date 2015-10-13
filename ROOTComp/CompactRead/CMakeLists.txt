if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-std=c++11 -pedantic-errors -Wall -Wextra -Wfatal-errors -Wwrite-strings -Woverloaded-virtual -fno-nonansi-builtins -fno-gnu-keywords -fstrict-aliasing")
endif()

include_directories(${PROJECT_SOURCE_DIR}/Support)
include_directories(${PROJECT_SOURCE_DIR}/../userinc)
add_library(CompactRead CompactIO.cpp CompactIO.h CompactImport.cpp CompactImport.h ROOTOutput.cpp ROOTOutput.h)