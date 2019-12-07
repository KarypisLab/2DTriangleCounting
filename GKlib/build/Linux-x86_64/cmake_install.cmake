# Install script for directory: /home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/karypisg/tomxx030/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/build/Linux-x86_64/libGKlib.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/GKlib.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/gk_arch.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/gk_defs.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/gk_externs.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/gk_getopt.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/gk_macros.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/gk_mkblas.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/gk_mkmemory.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/gk_mkpqueue.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/gk_mkpqueue2.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/gk_mkrandom.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/gk_mksort.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/gk_mkutils.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/gk_proto.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/gk_struct.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/gk_types.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/gkregex.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/ms_inttypes.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/ms_stat.h"
    "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/ms_stdint.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/build/Linux-x86_64/test/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/karypisg/tomxx030/distr_tc_cp/graphchallenge/distr_triangle_biggraphs/distr_triangle_2d/GKlib/build/Linux-x86_64/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
