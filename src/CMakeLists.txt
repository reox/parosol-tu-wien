cmake_minimum_required(VERSION 2.8)

PROJECT( parosol )
SET( parosol_src GReader.cpp GWriter.cpp JacobiSmoother.cpp Toolbox.cpp AsciiImage.cpp HDF5Image.cpp AsciiImageMirrored.cpp fem.cpp est_ev.cpp Chebyshev.cpp PCGSolver.cpp main.cpp)

ADD_EXECUTABLE( parosol ${parosol_src} )

#Debuggin
#SET( CMAKE_VERBOSE_MAKEFILE on )
set(CMAKE_CXX_FLAGS "-Wall ${CMAKE_CXX_FLAGS}")

# HDF5
# FIND_PACKAGE(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
set(LIBS ${LIBS} ${HDF5_LIBRARIES})

# Doxygen
include(FindDoxygen)
IF (DOXYGEN_FOUND)
  SET(DOXYGEN_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
  SET(DOXYGEN_OUTPUT_DIR doc)
  CONFIGURE_FILE(Doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  ADD_CUSTOM_TARGET(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
ELSE (DOXYGEN_FOUND)
    MESSAGE(STATUS "WARNING: Doxygen not found - Reference manual will not be created")
ENDIF (DOXYGEN_FOUND)

target_link_libraries(parosol ${LIBS})

install(TARGETS parosol DESTINATION bin)
