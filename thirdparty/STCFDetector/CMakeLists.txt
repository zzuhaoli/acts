cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

find_package ( DD4hep 1.21 REQUIRED COMPONENTS DDCore DDRec)
#include ( ${DD4hep_DIR}/cmake/DD4hep.cmake )
#include ( ${DD4hep_DIR}/cmake/DD4hepBuild.cmake )


if(TARGET XercesC::XercesC)
  SET(OPT_XERCESC XercesC::XercesC)
endif()

add_library(
  STCFDetector SHARED
  src/supportMDI_01.cpp 
  src/detectorCuRwell_02.cpp
  src/detectorMDC_03.cpp
  src/detectorRICH_01.cpp 
  src/detectorFTOF_01.cpp 
  src/detectorMUC_02.cpp 
  )

target_link_libraries(STCFDetector PUBLIC DD4hep::DDCore DD4hep::DDRec)
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})
set_target_properties(STCFDetector PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
dd4hep_set_version(STCFDetector MAJOR 1 MINOR 0 PATCH 0)
dd4hep_generate_rootmap(STCFDetector)


install(TARGETS STCFDetector LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})


configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/setup.sh.in
  ${PROJECT_BINARY_DIR}/this_stcf.sh
  @ONLY)
install(
  FILES ${PROJECT_BINARY_DIR}/this_stcf.sh
  DESTINATION ${CMAKE_INSTALL_BINDIR})

