cmake_minimum_required(VERSION 3.28)

add_executable(test test.cpp)

target_include_directories(
        test
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../../googletest/googlemock/include
)

target_include_directories(
        test
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../../googletest/googletest/include
)

target_include_directories(
        test
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

target_link_libraries(test gtest gtest_main)


set_target_properties(
        test PROPERTIES
        CXX_STANDARD 20
        CXX_STANDARD_REQUIRED ON
)