cmake_minimum_required(VERSION 3.1)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(utilities STATIC
    src/misc
)

SET_TARGET_PROPERTIES( utilities PROPERTIES COMPILE_FLAGS -O3)
