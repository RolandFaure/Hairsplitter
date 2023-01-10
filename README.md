# Hairsplitter

Splits contigs into their different haplotypes (or repeats into their different versions).

# What is Hairsplitter ?

`Hairsplitter` takes as input an assembly (obtained by any means) and the long reads (including high-error rate long reads) used to build this assembly. For each contig it checks if the contig was built using reads from different haplotypes/regions. If it was, `Hairsplitter` separates the reads into as many groups as necessary and computes the different versions (e.g. alleles) of the contig actually present in the genome. It outputs a new assembly, where different versions of contigs are not collapsed into one but assembled separately.

# Why is it useful ?

`Hairsplitter` is mainly useful if you are trying to obtain a phased assembly. The main advantage of `Hairsplitter` compared to other techniques is that it is totally parameters-free. Most importantly, it does not requires to know the ploidy of the organism, and can infer different ploidies corresponding to different contigs. It can thus be used just as well on haploid assemblies (to improve the assembly of duplications) as on complex allotetraploids (to assemble separately the haplotypes). Just run the assembly through!

# Installation

A conda package is in preparation but is not available yet. For now, it is necessary to download and compile the code.

## Dependancies

### Quick conda dependancies

You can create a conda environment with all dependencies installed by typing: 
```
conda create -c bioconda -c conda-forge -c anaconda -n hairsplitter minimap2 racon cmake gxx python=3.10.4 scipy numpy 
conda activate hairsplitter
```

### List of dependencies

- [minimap2](https://github.com/lh3/minimap2): this is not absolutely necessary. Another aligner can be used to manually align the reads.
- [racon](https://github.com/isovic/racon)
- CMake >= 3.8.12, make, gcc >= 11
- Python3 with numpy and scipy
- [GraphUnzip](https://github.com/nadegeguiglielmoni/GraphUnzip): this is included in the `Hairsplitter` folder in Hairsplitter/GraphUnzip/graphunzip.py, you do not need to install it separately. It is a python file, thus it will not work if taken out of its folder.
 
## Download & Compilation

To download and compile, run
```
git clone --recursive https://github.com/RolandFaure/Hairsplitter.git
cd Hairsplitter
mkdir build
cd build
cmake .. -DCMAKE_CXX_COMPILER=/path/to/g++ -DCMAKE_C_COMPILER=/path/to/gcc
make
```

Because of some optimization implemented directly in WFA2 (which we do not develop), Hairsplitter by default only runs on CPUs supporting AVX2 instructions (recent CPUs). To run HairSplitter on older CPUs, you must align the reads base-per-base to the assembly (e.g. using the -a option of minimap2) and then give the SAM file to HairSplitter through the -a option.

# Usage

## Quick start

Let's say you used `reads.fastq` (any long reads) to build the assembly `assembly.gfa` (with any assembler). To improve/phase your assembly using `Hairsplitter`, run
```
Hairsplitter -f reads.fastq -i assembly.gfa -o new_assembly.gfa
```

## Options

```bash
SYNOPSIS
        Hairsplitter -f <raw reads> -i <assembly> -o <output assembly> [-a
                                   <aligned reads>] [-q <output GAF>] [-p] [-t <threads>] [-s]
                                   [--path-to-minimap2 <path to minimap2>] [--path-to-miniasm <path
                                   to miniasm>] [--path-to-racon <path to racon>]
                                   [--path-to-graphunzip <path to graphunzip>]

OPTIONS
        -f, --fastq Sequencing reads
        -i, --assembly
                    Original assembly in GFA or FASTA format

        -o, --outputGFA
                    Output assembly file, same format as input

        -a, --aln-on-asm
                    Reads aligned on assembly (PAF or SAM format)

        -p, --polish
                    Use this option if the assembly is not polished

        -t, --threads
                    Number of threads

        -s, --dont_simplify
                    Don't rename the contigs and don't merge them

        --path-to-minimap2
                    Path to the executable minimap2 (if not in PATH)

        --path-to-racon
                    Path to the executable racon (if not in PATH)

        --path-to-graphunzip
                    Path to graphunzip.py (if not in PATH or in default folder)
                  
```

# Installation issues
 Most installation issues that we have seen yet stem from the use of too old compilers. g++ and gcc have to support c++17. Sometimes their default versions (especially on servers) are too old. Specify modern versions manually to cmake using `-DCMAKE_CXX_COMPILER=/path/to/modern/g++` and `-DCMAKE_C_COMPILER=/path/to/modern/gcc`.




