if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-std=c++11 -pedantic-errors -Wall -Wextra -Wfatal-errors -Wwrite-strings -Woverloaded-virtual -fno-nonansi-builtins -fno-gnu-keywords -fstrict-aliasing -g3")
endif()

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/CompactRead)
include_directories(${PROJECT_SOURCE_DIR}/Support)
add_library(Support OptionsParser.h OptionsParser.cpp ScanCuts.h ScanCuts.cpp)