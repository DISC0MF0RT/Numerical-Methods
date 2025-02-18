cmake_minimum_required(VERSION 3.10)

project(NM)

include_directories(
    ${PROJECT_SOURCE_DIR}/src
)

add_library(matrixLib STATIC
    src/matrix.cpp
)



set(CMAKE_CXX_STANDARD 17)

add_executable(NM src/main.cpp)

target_link_libraries(NM
    PRIVATE matrixLib
)
