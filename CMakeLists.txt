cmake_minimum_required(VERSION 2.8)
project(MPITC C)

set(GKLIB_PATH "GKlib" CACHE PATH "path to GKlib")
set(SHARED FALSE CACHE BOOL "build a shared library")

if(MSVC)
  set(ORION_INSTALL FALSE)
else()
  set(ORION_INSTALL TRUE)
endif()

# Add GKlib's cmake settings
include(${GKLIB_PATH}/GKlibSystem.cmake)

# Source files
file(GLOB mpitc_sources *.c)

# Include directories
include_directories(${GKLIB_PATH})
include_directories(.)
include_directories(~/local/include)

# Link directories 
link_directories(~/local/lib)

# Build
add_executable(mpitc ${GKlib_sources} ${mpitc_sources})

foreach(prog mpitc)
  target_link_libraries(${prog} m)
endforeach(prog)

# Install
if(ORION_INSTALL)
  install(TARGETS mpitc 
    RUNTIME DESTINATION bin)
endif()

# Try to find subversion revision.
set(SVNREV "")
file(TO_NATIVE_PATH ${PROJECT_SOURCE_DIR}/.svn svn_dir)
if(IS_DIRECTORY ${svn_dir})
  include(FindSubversion)
  if(Subversion_FOUND)
    Subversion_WC_INFO(${PROJECT_SOURCE_DIR} metis)
    set(SVNREV ${metis_WC_REVISION})
  endif(Subversion_FOUND)
endif()
add_definitions(-DSVNINFO="${SVNREV}")
