project(CRWAVES CXX)
cmake_minimum_required(VERSION 3.19)

set( CMAKE_VERBOSE_MAKEFILE off )
set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )
set( CMAKE_CXX_COMPILER "/usr/local/bin/g++-10" )
set( CMAKE_CXX_FLAGS "-O3 -fopenmp" )
set( CRWAVES_EXTRA_INCLUDES )
set( CRWAVES_EXTRA_LIBRARIES )

include_directories( "${CMAKE_SOURCE_DIR}" )

set( CRWAVES_SRCS 
     diffusion.cpp
     dump.cpp
     energy_losses.cpp
     evolutors.cpp
     evolve.cpp
     main.cpp
     params.cpp
     sources.cpp
     tridiag.cpp
     utilities.cpp
     waves.cpp
)

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
