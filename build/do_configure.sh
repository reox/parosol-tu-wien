#!/bin/bash
cmake \
    -DCMAKE_CXX_COMPILER:FILEPATH=mpic++ \
    -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo \
    -DHDF5_LIBRARIES:PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so \
    -DHDF5_INCLUDE_DIRS:PATH=/usr/include/hdf5/serial \
    -DCMAKE_CXX_FLAGS:STRING="-lz" \
    ../src
