set(NMCVT_TARGET "NMCVT")
file(GLOB_RECURSE NMCVT_HEADERS *.h)
file(GLOB_RECURSE NMCVT_SOURCES *.cpp)

add_library(${NMCVT_TARGET} SHARED ${NMCVT_SOURCES})
target_include_directories(${NMCVT_TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} Base)
target_link_libraries(${NMCVT_TARGET} igl::core geogram stan-math Base)