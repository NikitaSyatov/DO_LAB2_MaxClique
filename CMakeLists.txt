cmake_minimum_required(VERSION 3.10)
project(MaxCliqueSolver
    VERSION 1.0.0
    DESCRIPTION "C++ max clique solver"
    LANGUAGES CXX
)

# Настройка стандарта C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(BUILD_BIN_DIR ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BUILD_BIN_DIR})

# Опции компиляции
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-Wall -g -O0)
else()
    add_compile_options(-Wall -O3)
endif()

# Добавление исполняемого файла
add_executable(${PROJECT_NAME}
    source/clique_stud.cpp
)

# Установка свойств цели
set_target_properties(${PROJECT_NAME} PROPERTIES
    OUTPUT_NAME "MaxClique"
    VERSION ${PROJECT_VERSION}
)

# Установка путей для установки (опционально)
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/install/bin
)