cmake_minimum_required(VERSION 3.8.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(Hairsplitter)
#add_executable(${PROJECT_NAME} "main.cpp")

# Local header files here ONLY
set(TARGET_H
    split_reads.h
    split_reads_BiHap.h
    input_output.h
    modify_gfa.h
    robin_hood.h
    read.h
    sequence.h
    tools.h
    Partition.h
    clean_graph.h
    cluster_graph.h
    clipp.h
    call_variants.h
    phase_variants.h
    reassemble_unaligned_reads.h
    separate_reads.h
   )

# Local source files here
set(TARGET_SRC
    split_reads.cpp
    split_reads_BiHap.cpp
    edlib/src/edlib.cpp
    input_output.cpp
    modify_gfa.cpp
    read.cpp
    sequence.cpp
    tools.cpp
    Partition.cpp
    clean_graph.cpp
    cluster_graph.cpp
    call_variants.cpp
    phase_variants.cpp
    reassemble_unaligned_reads.cpp
    separate_reads.cpp
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

# file(GLOB SOURCES "main.cpp" "split_reads.cpp" "cluster_graph.cpp" "reassemble_unaligned_reads.cpp" "split_reads_BiHap.cpp" 
#     "input_output.cpp" "modify_gfa.cpp" "read.cpp" "sequence.cpp" "tools.cpp" "Partition.cpp" )
# add_executable(Hairsplitter ${SOURCES} edlib/src/edlib.cpp)
# # add_dependencies(Hairsplitter WFA)


# # target_link_libraries(Hairsplitter PRIVATE wfa2)#pthread)

# target_include_directories(Hairsplitter PRIVATE edlib/include)
# target_compile_options (Hairsplitter PRIVATE -g)
# target_compile_options (Hairsplitter PRIVATE -lm)
# target_compile_options (Hairsplitter PRIVATE -fopenmp)
# target_compile_options (Hairsplitter PRIVATE -std=c++17)
# target_compile_options (Hairsplitter PRIVATE -O3)
# target_compile_options (Hairsplitter PRIVATE -march=x86-64)

file(GLOB CONVERT_SOURCE "fa2gfa.cpp")
add_executable(fa2gfa ${CONVERT_SOURCE})

file (GLOB CONVERT_SOURCE2 "gfa2fa.cpp")
add_executable(gfa2fa ${CONVERT_SOURCE2})

file (GLOB SOURCE_CLEAN_GRAPH "clean_graph.cpp" "input_output.cpp" "tools.cpp")
add_executable(clean_graph ${SOURCE_CLEAN_GRAPH})
target_include_directories(clean_graph PRIVATE edlib/include)

file (GLOB SOURCE_CALL_VARIANTS "call_variants.cpp" "input_output.cpp" "robin_hood.h" "tools.cpp" "read.cpp" "sequence.cpp" "Partition.cpp")
add_executable(call_variants ${SOURCE_CALL_VARIANTS})
target_include_directories(call_variants PRIVATE edlib/include)
target_compile_options (call_variants PRIVATE -g)
target_compile_options (call_variants PRIVATE -lm)
target_compile_options (call_variants PRIVATE -fopenmp)
target_compile_options (call_variants PRIVATE -std=c++17)
target_compile_options (call_variants PRIVATE -O3)
target_compile_options (call_variants PRIVATE -march=x86-64)

file (GLOB SOURCE_FILTER_VARIANTS "filter_variants.cpp" "input_output.cpp" "robin_hood.h" "tools.cpp" "read.cpp" "sequence.cpp" "Partition.cpp")
add_executable(filter_variants ${SOURCE_FILTER_VARIANTS})
target_include_directories(filter_variants PRIVATE edlib/include)
target_compile_options (filter_variants PRIVATE -g)
target_compile_options (filter_variants PRIVATE -lm)
target_compile_options (filter_variants PRIVATE -fopenmp)
target_compile_options (filter_variants PRIVATE -std=c++17)
target_compile_options (filter_variants PRIVATE -O3)
target_compile_options (filter_variants PRIVATE -march=x86-64)

file (GLOB SOURCE_SEPARATE_READS "separate_reads.cpp" "input_output.cpp" "cluster_graph.cpp" "robin_hood.h" "tools.cpp" "read.cpp" "sequence.cpp" "Partition.cpp")
add_executable(separate_reads ${SOURCE_SEPARATE_READS})
target_include_directories(separate_reads PRIVATE edlib/include)
target_compile_options (separate_reads PRIVATE -g)
target_compile_options (separate_reads PRIVATE -lm)
target_compile_options (separate_reads PRIVATE -fopenmp)
target_compile_options (separate_reads PRIVATE -std=c++17)
target_compile_options (separate_reads PRIVATE -O3)
target_compile_options (separate_reads PRIVATE -march=x86-64)

file (GLOB SOURCE_CREATE_NEW_CONTIG "create_new_contigs.cpp" "input_output.cpp" "robin_hood.h" "tools.cpp" "read.cpp" "sequence.cpp" "Partition.cpp")
add_executable(create_new_contigs ${SOURCE_CREATE_NEW_CONTIG} edlib/src/edlib.cpp)
target_include_directories(create_new_contigs PRIVATE edlib/include)
target_compile_options (create_new_contigs PRIVATE -g)
target_compile_options (create_new_contigs PRIVATE -lm)
target_compile_options (create_new_contigs PRIVATE -fopenmp)
target_compile_options (create_new_contigs PRIVATE -std=c++17)
target_compile_options (create_new_contigs PRIVATE -O3)
target_compile_options (create_new_contigs PRIVATE -march=x86-64)


find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    # target_link_libraries(Hairsplitter PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(call_variants PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(filter_variants PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(separate_reads PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(create_new_contigs PRIVATE OpenMP::OpenMP_CXX)
endif()

#for OpenMP: https://answers.ros.org/question/64231/error-in-rosmake-rgbdslam_freiburg-undefined-reference-to-gomp/
