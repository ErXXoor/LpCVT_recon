cmake_minimum_required(VERSION 3.2)
project(LpCVT_recon)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(libigl)
igl_include(glfw)
include(geogram)

set(CMAKE_CXX_STANDARD 17)
add_subdirectory(src)
