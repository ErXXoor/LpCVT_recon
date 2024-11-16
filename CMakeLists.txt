cmake_minimum_required(VERSION 3.2)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(LpCVT_recon)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(stan_math)
include(libigl)
igl_include(glfw)
igl_include(copyleft comiso)
include(geogram)

find_package(PCL REQUIRED COMPONENTS common io filters)

set(CMAKE_CXX_STANDARD 11)
add_subdirectory(Base)

set(COMPILE_NMCVT_MESH OFF)
set(COMPILE_VIS OFF)
set(COMPILE_NMCVT_PC ON)

if (COMPILE_NMCVT_MESH)
    add_subdirectory(src)
endif ()
if (COMPILE_VIS)
    add_subdirectory(vis)
endif ()
if (COMPILE_NMCVT_PC)
    add_subdirectory(src_pc)
endif ()