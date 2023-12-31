# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-src"
  "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-build"
  "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-subbuild/glad-populate-prefix"
  "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-subbuild/glad-populate-prefix/tmp"
  "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-subbuild/glad-populate-prefix/src/glad-populate-stamp"
  "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-subbuild/glad-populate-prefix/src"
  "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-subbuild/glad-populate-prefix/src/glad-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-subbuild/glad-populate-prefix/src/glad-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-subbuild/glad-populate-prefix/src/glad-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
