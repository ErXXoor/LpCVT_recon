# Install script for directory: /Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/geogram-src/src/lib/geogram_gfx

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

if(CMAKE_INSTALL_COMPONENT STREQUAL "devkit" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/lib/libgeogram_gfx.1.8.6.dylib"
    "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/lib/libgeogram_gfx.1.dylib"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram_gfx.1.8.6.dylib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram_gfx.1.dylib"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      execute_process(COMMAND "/usr/bin/install_name_tool"
        -id "/usr/local/lib/libgeogram_gfx.1.dylib"
        -change "@rpath/libgeogram.1.dylib" "/usr/local/lib/libgeogram.1.dylib"
        "${file}")
      execute_process(COMMAND /usr/bin/install_name_tool
        -delete_rpath "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/geogram-src/cmake-build-release/lib"
        -delete_rpath "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/lib"
        -add_rpath "/usr/local/lib"
        "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "${file}")
      endif()
    endif()
  endforeach()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "devkit" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/lib/libgeogram_gfx.dylib")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "devkit-full" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/lib/libgeogram_gfx.1.8.6.dylib"
    "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/lib/libgeogram_gfx.1.dylib"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram_gfx.1.8.6.dylib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram_gfx.1.dylib"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      execute_process(COMMAND "/usr/bin/install_name_tool"
        -id "/usr/local/lib/libgeogram_gfx.1.dylib"
        -change "@rpath/libgeogram.1.dylib" "/usr/local/lib/libgeogram.1.dylib"
        "${file}")
      execute_process(COMMAND /usr/bin/install_name_tool
        -delete_rpath "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/geogram-src/cmake-build-release/lib"
        -delete_rpath "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/lib"
        -add_rpath "/usr/local/lib"
        "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "${file}")
      endif()
    endif()
  endforeach()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "devkit-full" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/lib/libgeogram_gfx.dylib")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "devkit" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/geogram1/geogram_gfx" TYPE DIRECTORY FILES "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/geogram-src/src/lib/geogram_gfx/." FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "devkit-full" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/geogram1/geogram_gfx" TYPE DIRECTORY FILES "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/geogram-src/src/lib/geogram_gfx/." FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/geogram-build/geogram_gfx1.pc")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/_deps/geogram-build/src/lib/geogram_gfx/third_party/cmake_install.cmake")

endif()

