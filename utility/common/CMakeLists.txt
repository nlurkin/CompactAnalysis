if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-g -std=c++11 -pedantic-errors -Wall -Wextra -Werror -Wfatal-errors -Wwrite-strings -Woverloaded-virtual -fno-nonansi-builtins -fno-gnu-keywords -fstrict-aliasing")
endif()

add_library(common pi0DalitzHeader.h style.cpp ConfigFile.h ConfigFile.cpp Fitter.h Fitter.cpp 
        RunWeights.cpp RunWeights.h Sample.cpp Sample.h FitMCSample.cpp FitMCSample.h FitDataSample.cpp
        FitDataSample.h)