cmake_minimum_required(VERSION 3.2)
project(LpCVT_recon)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(stan_math)
include(libigl)
igl_include(glfw)
igl_include(copyleft comiso)
include(geogram)

set(CMAKE_CXX_STANDARD 11)
add_subdirectory(Base)
add_subdirectory(src)
add_subdirectory(vis)