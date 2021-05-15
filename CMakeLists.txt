cmake_minimum_required(VERSION 3.17)
project(RT)

set(CMAKE_CXX_STANDARD 20)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O1")

set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--stack,999999")
add_subdirectory(extern/glfw-3.3.4)
find_package(OpenGL REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(include)
add_executable(RT src/main.cpp src/Ray.cpp include/Ray.h src/Hitable.cpp include/Hitable.h src/Camera.cpp include/Camera.h src/Material.cpp include/Material.h)
target_link_libraries(RT glfw)
target_link_libraries(RT OpenGL::GL)
target_link_libraries (RT Eigen3::Eigen)
