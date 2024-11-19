set(BASE_TARGET "Base")
file(GLOB_RECURSE BASE_HEADERS *.h)
file(GLOB_RECURSE BASE_SOURCES *.cpp)

add_library(${BASE_TARGET} SHARED ${BASE_SOURCES}
        Base/NNSearch.h
        NNSearch.cpp)
target_include_directories(${BASE_TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${BASE_TARGET} igl::glfw geogram igl_copyleft::comiso)