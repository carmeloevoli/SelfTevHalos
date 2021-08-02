project(CRWAVES CXX)
cmake_minimum_required(VERSION 3.1)

set( CMAKE_VERBOSE_MAKEFILE off )
set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )
set( CMAKE_CXX_COMPILER "/usr/local/bin/g++-10" )
set( CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -std=c++14 -fopenmp -DOPENMP -DTHREADS=16" )
set( CRWAVES_INCLUDES )
set( CRWAVES_LIBRARIES )

# include main directory
include_directories( "${CMAKE_SOURCE_DIR}" )

# create the output dir
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# utilities (provided)
add_subdirectory(libs/utilities)
list(APPEND CRWAVES_LIBRARIES utilities)
list(APPEND CRWAVES_INCLUDES libs/utilities/include)

# GSL (external)
find_package(GSL)
if(GSL_FOUND)
  set(HAVE_LIBGSL 1)
  list(APPEND CRWAVES_INCLUDES ${GSL_INCLUDE_DIR})
  list(APPEND CRWAVES_LIBRARIES ${GSL_LIBRARIES})
else(GSL_FOUND)
  message(FATAL_ERROR "Require GSL. Set GSL_HOME")
endif(GSL_FOUND)

# add source files
set( CRWAVES_SRCS 
    src/builders.cpp
    src/common.cpp
    src/diffusion.cpp
    src/dump.cpp
    src/evolutors.cpp
    src/evolve.cpp
    src/params.cpp
    src/tridiag.cpp
    src/waves.cpp
)

# make library
include_directories (include ${CRWAVES_INCLUDES})
add_library( CRWAVES_LIB ${CRWAVES_SRCS} ${CRWAVES_INCLUDES} )

# MAKE EXECUTABLES
add_executable(crwaves src/main.cpp)
target_link_libraries (crwaves CRWAVES_LIB ${CRWAVES_LIBRARIES})

add_executable(test_common tests/test_common.cpp)
target_link_libraries (test_common CRWAVES_LIB ${CRWAVES_LIBRARIES})

add_executable(test_nintegrate tests/test_nintegrate.cpp)
target_link_libraries (test_nintegrate CRWAVES_LIB ${CRWAVES_LIBRARIES})


# Show summary.
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "CRWAVES_INCLUDES ........ = ${CRWAVES_INCLUDES}")
message (STATUS "CRWAVES_LIBRARIES ....... = ${CRWAVES_LIBRARIES}")