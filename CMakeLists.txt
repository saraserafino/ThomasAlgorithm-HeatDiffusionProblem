cmake_minimum_required(VERSION 3.12)

project(project VERSION 1.0
    DESCRIPTION "exam_project_Serafino"
    )

add_compile_options(-std=c++17)

# Search for muparserx
find_package(muparserx 4.0 REQUIRED)
include_directories(${muparserx_INCLUDE_DIRS})

# Search for Eigen
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})

# Search for Pybind
find_package(pybind11 REQUIRED)
include_directories(SYSTEM ${pybind11_INCLUDE_DIRS})

# Add the Python modules
pybind11_add_module(moduleH
    MatrixModule/src/Matrix.cpp
    MatrixModule/include/ThomasSolver.tpl.hpp
    MatrixModule/src/HeatDiffusion.cpp
    MatrixModule/src/Matrix_py.cpp
    )

target_link_libraries(moduleH PRIVATE ${MUPARSERX_LIBRARY})
target_include_directories(moduleH PRIVATE MatrixModule/include)