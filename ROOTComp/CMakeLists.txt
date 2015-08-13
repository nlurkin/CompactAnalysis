cmake_minimum_required (VERSION 2.6 FATAL_ERROR)

project (pi0d_selection)

#require boost
set(BOOST_ROOT $ENV{BOOST})
find_package(Boost 1.53.0 REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})
  
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

# Own libraries
add_subdirectory(CompactRead)
add_subdirectory(Support)

include_directories(CompactRead)
include_directories(Support)

add_executable(pi0d_selection pi0d_Selection.cpp)
add_executable(alt_id pid.cpp)
add_executable(alt_id_kmu3 pid_mu.cpp)

target_link_libraries(pi0d_selection CompactRead Support) #Own
target_link_libraries(pi0d_selection exportClasses funLib) #Compact
target_link_libraries(pi0d_selection ${Boost_LIBRARIES} ${ROOT_LIBRARIES}) #External

target_link_libraries(alt_id CompactRead Support) #Own
target_link_libraries(alt_id exportClasses funLib) #Compact
target_link_libraries(alt_id ${Boost_LIBRARIES} ${ROOT_LIBRARIES}) #External

target_link_libraries(alt_id_kmu3 CompactRead Support) #Own
target_link_libraries(alt_id_kmu3 exportClasses funLib) #Compact
target_link_libraries(alt_id_kmu3 ${Boost_LIBRARIES} ${ROOT_LIBRARIES}) #External

install(TARGETS pi0d_selection alt_id alt_id_kmu3 DESTINATION .)
