cmake_minimum_required(VERSION 2.8)

project(IMU)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_FLAGS "-std=c++11")

find_package(Eigen REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include
		    ${EIGEN_INCLUDE_DIRS})
add_subdirectory(src)
