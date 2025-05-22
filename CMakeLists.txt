cmake_minimum_required(VERSION 3.19)
PROJECT(Parosol VERSION 1.0.0)

SET(PAROSOL_BIN parosol)

# Usually this is what you want...
set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build.")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

## Library Dependencies
SET(HDF5_PREFER_PARALLEL TRUE)

FIND_PACKAGE(MPI REQUIRED)
FIND_PACKAGE(HDF5 REQUIRED)
FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)

# Check if we have indeed found the parallel version
if(NOT HDF5_IS_PARALLEL)
    MESSAGE(FATAL_ERROR "Parallel HDF5 library was not found! Have you installed the MPI variant of HDF5?")
endif()

# Binaries to create
ADD_EXECUTABLE(${PAROSOL_BIN}
    src/main.cpp
    src/GReader.cpp
    src/GWriter.cpp
    src/JacobiSmoother.cpp
    src/Toolbox.cpp
    src/AsciiImage.cpp
    src/HDF5Image.cpp
    src/AsciiImageMirrored.cpp
    src/fem.cpp
    src/est_ev.cpp
    src/Chebyshev.cpp
    src/PCGSolver.cpp
)

# add_definitions(-DSEQUENTIAL_HDF5)

# Available since cmake 3.9
target_link_libraries(${PAROSOL_BIN} PUBLIC MPI::MPI_CXX)
# Target for HDF5 is supported with cmake >= 3.19
target_link_libraries(${PAROSOL_BIN} PUBLIC HDF5::HDF5)
# Supported with Eigen 3.3
target_link_libraries(${PAROSOL_BIN} PUBLIC Eigen3::Eigen)

