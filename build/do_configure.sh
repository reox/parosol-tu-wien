#!/bin/bash
cmake \
    -DCMAKE_CXX_COMPILER:FILEPATH=mpic++ \
    -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo \
    ../src
