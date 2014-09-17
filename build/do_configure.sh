#!/bin/bash
/usr2/elan/installs/bin/cmake \
    -DCMAKE_CXX_COMPILER:FILEPATH=mpic++ \
    -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo \
    -DHDF5_LIBRARIES:PATH=/usr2/elan/installs2/lib/libhdf5.so \
    -DHDF5_INCLUDE_DIRS:PATH=/usr2/elan/installs2/include \
    -DCMAKE_CXX_FLAGS:STRING=-I/usr2/elan/parosol \
    ../src
