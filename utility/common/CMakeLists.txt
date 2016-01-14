if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-g -std=c++11 -pedantic-errors -Wall -Wextra -Werror -Wfatal-errors -Wwrite-strings -Woverloaded-virtual -fno-nonansi-builtins -fno-gnu-keywords -fstrict-aliasing")
endif()

set(drawerHdr Interface/HistoDrawer.h Drawer/InputFitDrawer.h Drawer/FitResultDrawer.h Drawer/CombineDrawer.h Drawer/Drawer.h Drawer/StackDrawer.h Drawer/StackRatioDrawer.h) 
set(drawerSrc Interface/HistoDrawer.cpp Drawer/InputFitDrawer.cpp Drawer/FitResultDrawer.cpp Drawer/CombineDrawer.cpp Drawer/Drawer.cpp Drawer/StackDrawer.cpp Drawer/StackRatioDrawer.cpp)

set(fitterHdr Fitter/Fitter.h Fitter/MinuitFitter.h Fitter/MinuitFitterNewROOT.h Fitter/MinuitFitterNewChi2.h)
set(fitterSrc Fitter/Fitter.cpp Fitter/MinuitFitter.cpp Fitter/MinuitFitterNewROOT.cpp Fitter/MinuitFitterNewChi2.cpp)

set(sampleHdr Interface/Sample.h Interface/DataSample.h Interface/MCSample.h Samples/FitMCSample.h Samples/FitDataSample.h)
set(sampleSrc Interface/Sample.cpp Interface/DataSample.cpp Interface/MCSample.cpp Samples/FitMCSample.cpp Samples/FitDataSample.cpp)

set(sampleHdr ${sampleHdr} Samples/CombineSample.h Samples/CombineMCSample.h Samples/CombineDataSample.h)
set(sampleSrc ${sampleSrc} Samples/CombineSample.cpp Samples/CombineMCSample.cpp Samples/CombineDataSample.cpp)

set(sampleHdr ${sampleHdr} Interface/SubSample.h)
set(sampleSrc ${sampleSrc} Interface/SubSample.cpp)

add_library(Drawer ${drawerHdr} ${drawerSrc})
add_library(Fitter ${fitterHdr} ${fitterSrc})
add_library(Combiner Combiner.h Combiner.cpp)
add_library(Samples ${sampleHdr} ${sampleSrc})

add_library(common pi0DalitzHeader.h style.cpp ConfigFile.h ConfigFile.cpp RunWeights.cpp RunWeights.h DataGetter.h DataGetter.cpp)