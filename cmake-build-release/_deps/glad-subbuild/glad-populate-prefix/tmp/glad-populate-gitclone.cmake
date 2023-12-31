# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if(EXISTS "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-subbuild/glad-populate-prefix/src/glad-populate-stamp/glad-populate-gitclone-lastrun.txt" AND EXISTS "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-subbuild/glad-populate-prefix/src/glad-populate-stamp/glad-populate-gitinfo.txt" AND
  "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-subbuild/glad-populate-prefix/src/glad-populate-stamp/glad-populate-gitclone-lastrun.txt" IS_NEWER_THAN "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-subbuild/glad-populate-prefix/src/glad-populate-stamp/glad-populate-gitinfo.txt")
  message(STATUS
    "Avoiding repeated git clone, stamp file is up to date: "
    "'/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-subbuild/glad-populate-prefix/src/glad-populate-stamp/glad-populate-gitclone-lastrun.txt'"
  )
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-src"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-src'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/local/bin/git"
            clone --no-checkout --config "advice.detachedHead=false" "https://github.com/libigl/libigl-glad.git" "glad-src"
    WORKING_DIRECTORY "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps"
    RESULT_VARIABLE error_code
  )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once: ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/libigl/libigl-glad.git'")
endif()

execute_process(
  COMMAND "/usr/local/bin/git"
          checkout "ead2d21fd1d9f566d8f9a9ce99ddf85829258c7a" --
  WORKING_DIRECTORY "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-src"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'ead2d21fd1d9f566d8f9a9ce99ddf85829258c7a'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/usr/local/bin/git" 
            submodule update --recursive --init 
    WORKING_DIRECTORY "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-src"
    RESULT_VARIABLE error_code
  )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-src'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-subbuild/glad-populate-prefix/src/glad-populate-stamp/glad-populate-gitinfo.txt" "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-subbuild/glad-populate-prefix/src/glad-populate-stamp/glad-populate-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/glad-subbuild/glad-populate-prefix/src/glad-populate-stamp/glad-populate-gitclone-lastrun.txt'")
endif()
