cmake_minimum_required(VERSION 3.8.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(GenomeTailor)

# Local header files here ONLY
set(TARGET_H
    scaffold.h,
    align.h,
    robin_hood.h,
    clipp.h
   )

# Local source files here
set(TARGET_SRC
    scaffold.cpp,
    align.cpp
    )

file(GLOB SOURCE "scaffold.cpp" "align.cpp")
add_executable(HS_GenomeTailor ${SOURCE})
target_compile_options (HS_GenomeTailor PRIVATE -g)
target_compile_options (HS_GenomeTailor PRIVATE -lm)
target_compile_options (HS_GenomeTailor PRIVATE -std=c++17)
target_compile_options (HS_GenomeTailor PRIVATE -O3)
target_compile_options (HS_GenomeTailor PRIVATE -march=x86-64)
