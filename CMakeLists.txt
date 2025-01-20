cmake_minimum_required(VERSION 3.2)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(LpCVT_recon)

set(CMAKE_CXX_STANDARD 11)
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(stan_math)
include(libigl)
igl_include(glfw)
include(geogram)

add_subdirectory(Base)
add_subdirectory(src)
#add_subdirectory(vis)