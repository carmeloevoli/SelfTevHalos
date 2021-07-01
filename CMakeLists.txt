project(CRWAVES CXX)
cmake_minimum_required(VERSION 3.1)

set( CMAKE_VERBOSE_MAKEFILE off )
set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )
set( CMAKE_CXX_COMPILER "/usr/local/bin/g++-10" )
set( CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -std=c++11 -fopenmp -DTHREADS=16" )
set( CRWAVES_EXTRA_INCLUDES )
set( CRWAVES_EXTRA_LIBRARIES )

include_directories( "${CMAKE_SOURCE_DIR}" )

set( CRWAVES_SRCS 
     src/diffusion.cpp
     src/dump.cpp
     src/energy_losses.cpp
     src/evolutors.cpp
     src/evolve.cpp
     src/main.cpp
     src/params.cpp
     src/sources.cpp
     src/tridiag.cpp
     src/utilities.cpp
     src/waves.cpp
)

# create the output dir
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# GSL (external)
find_package(GSL)
if(GSL_FOUND)
  set(HAVE_LIBGSL 1)
  list(APPEND CRWAVES_EXTRA_INCLUDES ${GSL_INCLUDE_DIR})
  list(APPEND CRWAVES_EXTRA_LIBRARIES ${GSL_LIBRARIES})
else(GSL_FOUND)
  message(FATAL_ERROR "Require GSL. Set GSL_HOME")
endif(GSL_FOUND)

# MAKE EXECUTABLE
include_directories (include ${CRWAVES_EXTRA_INCLUDES})

add_executable(CRWAVES ${CRWAVES_SRCS} )
target_link_libraries (CRWAVES ${CRWAVES_EXTRA_LIBRARIES})

# Show summary.
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "CRWAVES_INCLUDES ........ = ${CRWAVES_EXTRA_INCLUDES}")
message (STATUS "CRWAVES_LIBRARIES ....... = ${CRWAVES_EXTRA_LIBRARIES}")