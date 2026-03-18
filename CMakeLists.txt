cmake_minimum_required(VERSION 3.16)
project(nosehoover-chain LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Найдём библиотеки
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(GSL REQUIRED)

# Собираем всё в один бинарник
add_executable(nosehoover-chain
    eam.cpp
    export_xyz.cpp
    potential.cpp
    verlet_solve.cpp

    consts.h
    eam.h
    export_xyz.h
    potential.h
    verlet_solve.h
)

# Подключаем заголовки (все лежат рядом)
target_include_directories(nosehoover-chain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Линкуем зависимости
target_link_libraries(nosehoover-chain
    PRIVATE
        Eigen3::Eigen
        GSL::gsl
        GSL::gslcblas
)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3 -march=native -ffast-math -fno-math-errno -DNDEBUG)
endif()


target_compile_options(nosehoover-chain PRIVATE
    -Wall -Wextra -Wpedantic
    $<$<CONFIG:Release>:-O3 -march=native -ffast-math -fno-math-errno>
)
