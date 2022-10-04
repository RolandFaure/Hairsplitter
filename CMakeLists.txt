cmake_minimum_required(VERSION 3.8.12)

project(Hairsplitter)
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
add_executable(Hairsplitter ${SOURCES} edlib/src/edlib.cpp)
add_dependencies(Hairsplitter WFA)

target_include_directories(Hairsplitter PRIVATE edlib/include)
target_compile_options (Hairsplitter PRIVATE -fopenmp)
target_compile_options (Hairsplitter PRIVATE -std=c++17)
target_compile_options (Hairsplitter PRIVATE -O3)
target_compile_options (Hairsplitter PRIVATE -no-pie)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(Hairsplitter PRIVATE OpenMP::OpenMP_CXX)
endif()
target_link_libraries(Hairsplitter PRIVATE wfa2 pthread)

#for OpenMP: https://answers.ros.org/question/64231/error-in-rosmake-rgbdslam_freiburg-undefined-reference-to-gomp/
#set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS}")
#set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /DELAYLOAD:${OPENMP_DLL}.dll")

# #check for OpenMP
# find_package(OpenMP)
# if(OPENMP_FOUND)
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#   if(MSVC90 OR MSVC10)
#     if(MSVC90)
#       set(OPENMP_DLL VCOMP90)
#     elseif(MSVC10)
#       set(OPENMP_DLL VCOMP100)
#     endif(MSVC90)
#     set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} /DELAYLOAD:${OPENMP_DLL}D.dll")
#     set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /DELAYLOAD:${OPENMP_DLL}.dll")
#   endif()
# endif()



#if (NOT TARGET spoa)
#  add_subdirectory(spoa EXCLUDE_FROM_ALL)
#endif ()
#target_link_libraries(Hairsplitter spoa::spoa)
