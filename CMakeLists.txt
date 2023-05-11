cmake_minimum_required(VERSION 3.8.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(Hairsplitter)
#add_executable(${PROJECT_NAME} "main.cpp")

# Local header files here ONLY
set(TARGET_H
    split_reads.h
    input_output.h
    modify_gfa.h
    robin_hood.h
    read.h
    sequence.h
    tools.h
    Partition.h
    cluster_graph.h
    clipp.h
    phase_variants.h
    reassemble_unaligned_reads.h
   )

# Local source files here
set(TARGET_SRC
    split_reads.cpp
    edlib/src/edlib.cpp
    input_output.cpp
    modify_gfa.cpp
    read.cpp
    sequence.cpp
    tools.cpp
    Partition.cpp
    cluster_graph.cpp
    phase_variants.cpp
    reassemble_unaligned_reads.cpp
    )

# include(ExternalProject)


# ExternalProject_Add(
#     WFA
#     PREFIX WFA-ext
#     SOURCE_DIR ${CMAKE_SOURCE_DIR}/WFA2-lib/build
#     BUILD_IN_SOURCE true
#     LOG_BUILD ON
#     CONFIGURE_COMMAND ""
#     BUILD_COMMAND make
#     INSTALL_COMMAND ""
# )

# add_library(wfa2 INTERFACE)
# target_include_directories(wfa2 INTERFACE ${CMAKE_SOURCE_DIR}/WFA2-lib)
# target_link_directories(wfa2 INTERFACE ${CMAKE_SOURCE_DIR}/WFA2-lib/lib)
# target_link_libraries(wfa2 INTERFACE wfa wfacpp rt)

file(GLOB SOURCES "*.cpp")
add_executable(Hairsplitter ${SOURCES} edlib/src/edlib.cpp)
# add_dependencies(Hairsplitter WFA)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(Hairsplitter PRIVATE OpenMP::OpenMP_CXX)
endif()

# target_link_libraries(Hairsplitter PRIVATE wfa2)#pthread)

target_include_directories(Hairsplitter PRIVATE edlib/include)
target_compile_options (Hairsplitter PRIVATE -g)
target_compile_options (Hairsplitter PRIVATE -lm)
target_compile_options (Hairsplitter PRIVATE -fopenmp)
target_compile_options (Hairsplitter PRIVATE -std=c++17)
target_compile_options (Hairsplitter PRIVATE -O3)
target_compile_options (Hairsplitter PRIVATE -march=x86-64)

#for OpenMP: https://answers.ros.org/question/64231/error-in-rosmake-rgbdslam_freiburg-undefined-reference-to-gomp/

