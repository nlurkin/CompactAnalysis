if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-std=c++11 -pedantic-errors -Wall -Wextra -Wfatal-errors -Wwrite-strings -Woverloaded-virtual -fno-nonansi-builtins -fno-gnu-keywords -fstrict-aliasing -g3")
endif()

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/CompactRead)
include_directories(${PROJECT_SOURCE_DIR}/Support)

ROOT_GENERATE_DICTIONARY(ScanCuts.h LinkDef.h ScanCuts_DICT.cc "${ROOT_INCLUDE_DIR};${CMAKE_CURRENT_SOURCE_DIR};${PROJECT_SOURCE_DIR}/../userinc")
add_library(Support OptionsParser.h OptionsParser.cpp ScanCuts.h ScanCuts.cpp pid_res.h pid_res.cpp ScanCuts_DICT.cc)
add_library(Cuts SHARED ScanCuts.h ScanCuts.cpp ScanCuts_DICT.cc)
install(TARGETS Support Cuts DESTINATION lib)
install(FILES ScanCuts.h DESTINATION includes)
