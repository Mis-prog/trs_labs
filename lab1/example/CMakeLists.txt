cmake_minimum_required(VERSION 3.28)

add_executable(lab1_example main.cpp)

set_target_properties(
        lab1_example PROPERTIES
        CXX_STANDARD 20
        CXX_STANDARD_REQUIRED ON
)

target_include_directories(
        lab1_example
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
)