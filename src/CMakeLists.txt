# CMakeLists.txt for DataTree package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 2.4 FATAL_ERROR)
project(DataTree)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

ROOT_GENERATE_DICTIONARY(G__DataTreeEvent DataTreeConverterNA61.h DataTreeEvent.h LINKDEF LinkDef/DataTreeEventLinkDef.h)

#---Create a shared library with geneated dictionary
#add_library(DataTreeTrack SHARED DataTreeTrack.cxx)
add_library(DataTree SHARED 
DataTreeEvent.cxx
DataTreeTrack.cxx
DataTreePSDSection.cxx
DataTreePSDModule.cxx
DataTreeTOFHit.cxx
DataTreeV0Candidate.cxx
DataTreeMCTrack.cxx
DataTreeTrigger.cxx
DataTreeBPD.cxx
DataTreeConverterNA61.cxx
G__DataTreeEvent.cxx)

target_link_libraries(DataTree ${ROOT_LIBRARIES})

#---Create  a main program using the library
#add_executable(Main MainEvent.cxx)
#target_link_libraries(Main Event)
