# Install script for directory: /Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/geogram-src/src/tests

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
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

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/geogram-build/src/tests/test_nn_search/cmake_install.cmake")
  include("/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/geogram-build/src/tests/test_convex_cell/cmake_install.cmake")
  include("/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/geogram-build/src/tests/bench_load/cmake_install.cmake")
  include("/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/geogram-build/src/tests/test_locks/cmake_install.cmake")
  include("/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/geogram-build/src/tests/test_expansion_nt/cmake_install.cmake")
  include("/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/geogram-build/src/tests/test_HLBFGS/cmake_install.cmake")
  include("/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/geogram-build/src/tests/test_RVC/cmake_install.cmake")
  include("/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/geogram-build/src/tests/test_logger/cmake_install.cmake")
  include("/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/geogram-build/src/tests/test_CDT_2d/cmake_install.cmake")

endif()

