project(vis)
set(CMAKE_CXX_STANDARD 11)
file(GLOB_RECURSE SOURCE_FILES *.c *.cpp)
file(GLOB_RECURSE HEADER_FILES *.h *.hpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES} ${HEADER_FILES})

target_include_directories(${PROJECT_NAME} PRIVATE ${HEADER_FILES} Base)
target_link_libraries(${PROJECT_NAME} Base)