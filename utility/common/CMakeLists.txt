if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-g -std=c++11 -pedantic-errors -Wall -Wextra -Werror -Wfatal-errors -Wwrite-strings -Woverloaded-virtual -fno-nonansi-builtins -fno-gnu-keywords -fstrict-aliasing")
endif()

add_library(Drawer Drawer/HistoDrawer.h Drawer/HistoDrawer.cpp Drawer/InputFitDrawer.h Drawer/InputFitDrawer.cpp Drawer/FitResultDrawer.h Drawer/FitResultDrawer.cpp)
add_library(Fitter Fitter/Fitter.h Fitter/Fitter.cpp Fitter/MinuitFitter.h Fitter/MinuitFitter.cpp Fitter/MinuitFitterNewROOT.h Fitter/MinuitFitterNewROOT.cpp Fitter/MinuitFitterNewChi2.h Fitter/MinuitFitterNewChi2.cpp)
add_library(Samples Samples/Sample.h Samples/Sample.cpp Samples/FitMCSample.h Samples/FitMCSample.cpp Samples/FitDataSample.h Samples/FitDataSample.cpp)

add_library(common pi0DalitzHeader.h style.cpp ConfigFile.h ConfigFile.cpp RunWeights.cpp RunWeights.h)