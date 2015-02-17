if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-std=c++11 -pedantic-errors -Wall -Wextra -Werror -Wfatal-errors -Wwrite-strings -Woverloaded-virtual -fno-nonansi-builtins -fno-gnu-keywords -fstrict-aliasing")
endif()

add_library(CompactRead CompactIO.cpp CompactIO.h CompactImport.cpp CompactImport.h ROOTOutput.cpp ROOTOutput.h)