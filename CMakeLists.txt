project(RDKitjs)
cmake_minimum_required(VERSION 3.0)

set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.50.0 COMPONENTS ) 


set(CMAKE_C_FLAGS "-O2")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS}")

add_subdirectory(src)
