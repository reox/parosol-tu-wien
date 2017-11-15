cmake_minimum_required(VERSION 3.0)

add_subdirectory(src)
add_subdirectory(tools)

if (ENABLE_DOC)
    add_subdirectory(doc)
endif (ENABLE_DOC)
