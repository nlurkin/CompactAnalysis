cmake_minimum_required (VERSION 2.6 FATAL_ERROR)

project (pi0d_tools)

#require ROOT
include(${PROJECT_SOURCE_DIR}/FindROOT.cmake)
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

if(CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_CXX_FLAGS "-std=c++11 -pedantic-errors -Wall -Wextra -Wwrite-strings -Woverloaded-virtual -fno-nonansi-builtins -fno-gnu-keywords -fstrict-aliasing -g3")
endif()

#Compact
include_directories(${PROJECT_SOURCE_DIR}/../userinc)
link_directories(${PROJECT_SOURCE_DIR}/../obj)

include_directories(common)
add_subdirectory(common)

add_executable(combine batch/combine_batch.cpp)
add_executable(fit batch/fit_batch.cpp)
add_executable(make_bins batch/make_bins.cpp)
add_executable(prepareWeights batch/prepareWeights.cpp)
add_executable(toyMC batch/toyMC.cpp)

target_link_libraries(combine common exportClasses ${ROOT_LIBRARIES})
target_link_libraries(fit common exportClasses ${ROOT_LIBRARIES} Minuit)
target_link_libraries(make_bins common exportClasses ${ROOT_LIBRARIES})
target_link_libraries(prepareWeights common exportClasses ${ROOT_LIBRARIES})
target_link_libraries(toyMC common exportClasses ${ROOT_LIBRARIES} Minuit) 

install(TARGETS combine DESTINATION .)
install(TARGETS fit DESTINATION .)
install(TARGETS make_bins DESTINATION .)
install(TARGETS prepareWeights DESTINATION .)
install(TARGETS toyMC DESTINATION .)