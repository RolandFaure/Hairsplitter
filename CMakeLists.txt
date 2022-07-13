cmake_minimum_required(VERSION 3.8.12)

project(OverlapCheck)
#add_executable(${PROJECT_NAME} "main.cpp")

# Local header files here ONLY
set(TARGET_H
    check_overlaps.h
    input_output.h
    modify_gfa.h
    robin_hood.h
    read.h
    sequence.h
    tools.h
    Partition.h
    cluster_graph.h
    clipp.h
    Variant.h
    phase_variants.h
   )

# Local source files here
set(TARGET_SRC
    check_overlaps.cpp
    edlib/src/edlib.cpp
    input_output.cpp
    modify_gfa.cpp
    read.cpp
    sequence.cpp
    tools.cpp
    Partition.cpp
    cluster_graph.cpp
    Variant.cpp
    phase_variants.cpp
    )

include(ExternalProject)

ExternalProject_Add(
    WFA
    PREFIX WFA-ext
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/WFA2-lib
    BUILD_IN_SOURCE true
    LOG_BUILD ON
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make
    INSTALL_COMMAND ""
)

add_library(wfa2 INTERFACE)
target_include_directories(wfa2 INTERFACE ${CMAKE_SOURCE_DIR}/WFA2-lib)
target_link_directories(wfa2 INTERFACE ${CMAKE_SOURCE_DIR}/WFA2-lib/lib)
target_link_libraries(wfa2 INTERFACE wfa wfacpp)



file(GLOB SOURCES "*.cpp")
add_executable(OverlapCheck ${SOURCES} edlib/src/edlib.cpp)
add_dependencies(OverlapCheck WFA)

target_include_directories(OverlapCheck PRIVATE edlib/include)
target_compile_options (OverlapCheck PRIVATE -std=c++17)
target_compile_options (OverlapCheck PRIVATE -O3)
target_link_libraries(OverlapCheck wfa2 pthread)

#if (NOT TARGET spoa)
#  add_subdirectory(spoa EXCLUDE_FROM_ALL)
#endif ()
#target_link_libraries(OverlapCheck spoa::spoa)
