cmake_minimum_required(VERSION 2.8.12)

project(KMC)


# Set some default compile flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )


add_subdirectory(kmc_api)
