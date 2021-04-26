cmake_minimum_required(VERSION 3.15)
project(284Final)

find_package(glfw3 REQUIRED)
find_package(OpenGL REQUIRED)


set(CMAKE_CXX_STANDARD 11)

add_executable(284Final
        src/main.cpp
        src/Grid.cpp
        src/Matrix2D.cpp
        src/Particle.cpp
        src/Scene.cpp
        src/Simulator.cpp
        src/Snow.cpp
        src/Vector2D.cpp)

#add_library(A SHARED IMPORTED)
#set_property(TARGET A PROPERTY IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/GLFW/lib-mingw-w64/glfw3.dll)
target_link_libraries(284Final OpenGL::GL)
target_link_libraries(284Final glfw)