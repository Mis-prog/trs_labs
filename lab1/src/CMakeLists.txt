cmake_minimum_required(VERSION 3.10)

file(GLOB SRCS *.cpp)

add_library(
        lab1
        STATIC
        ${SRCS}
)

target_include_directories(
        lab1
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

set_target_properties(
        lab1 PROPERTIES
        CXX_STANDARD 20
        CXX_STANDARD_REQUIRED ON
)