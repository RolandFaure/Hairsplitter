# Hairsplitter

Splits contigs into their different haplotypes (or repeats into their different versions).

# What is Hairsplitter ?

`Hairsplitter` takes as input an assembly (obtained by any means) and the long reads (including high-error rate long reads) used to build this assembly. For each contig it checks if the contig was built using reads from different haplotypes/regions. If it was, `Hairsplitter` separates the reads into as many groups as necessary and computes the different versions (e.g. alleles) of the contig actually present in the genome. It outputs a new assembly, where different versions of contigs are not collapsed into one but assembled separately.

# Why is it useful ?

`Hairsplitter` can be used to refine a metagenomic assembly. Assemblers commonly collapse closely related strains as on single genome. HairSplitter can recover the lost strains.
`HairSplitter` is also useful for single-organism assembly, especially if you are trying to obtain a phased assembly. The main advantage of `Hairsplitter` compared to other techniques is that it is totally parameters-free. Most importantly, it does not requires to know the ploidy of the organism, and can infer different ploidies corresponding to different contigs. It can thus be used just as well on haploid assemblies (to improve the assembly of duplications) as on complex allotetraploids (to assemble separately the haplotypes). Just run the assembly through!

# Installation

A conda package is in preparation but is not available yet. For now, it is necessary to download and compile the code.

## Dependancies

### Quick conda dependancies

You can create and activate a conda environment with all dependencies installed by typing: 
```
conda create -c bioconda -c conda-forge -c anaconda -n hairsplitter minimap2 racon cmake gxx gcc python scipy numpy 
conda activate hairsplitter
```

### List of dependencies

- [minimap2](https://github.com/lh3/minimap2)
- [racon](https://github.com/isovic/racon)
- CMake >= 3.8.12, make, gcc >= 11, g++ >= 11
- Python3 with numpy and scipy

If Minimap2 and Racon are not in the PATH, their location should be specified through the `--path-to-minimap2` and `--path-to-racon` options.
 
## Download & Compilation

To download and compile, run
```
git clone https://github.com/RolandFaure/Hairsplitter.git
cd Hairsplitter
mkdir build
cd build
cmake .. -DCMAKE_CXX_COMPILER=/path/to/g++ -DCMAKE_C_COMPILER=/path/to/gcc
make
```

# Usage

## Quick start

Let's say you used `reads.fastq` (any long reads) to build the assembly `assembly.gfa` (with any assembler)(the assembly can be in gfa or fasta format). To improve/phase your assembly using `Hairsplitter`, run
```
Hairsplitter -f reads.fastq -i assembly.gfa -o hairsplitter_out/
```

In the folder hairsplitter_out, you will find the new assembly, named `hairsplitter\_assembly.gfa`. Another generated file is `hairsplitter\_summary.txt`, in which are written which contigs are duplicated and merged.

## Options

```bash
SYNOPSIS
        Hairsplitter -f <raw reads> -i <assembly> -o <output assembly> [-a
                                   <aligned reads>] [-q <output GAF>] [-p] [-t <threads>] [-s]
                                   [--path-to-minimap2 <path to minimap2>] [--path-to-miniasm <path
                                   to miniasm>] [--path-to-racon <path to racon>]
                                   [--path-to-graphunzip <path to graphunzip>]

OPTIONS
        -f, --fastq Sequencing reads (required)
        -i, --assembly
                    Original assembly in GFA or FASTA format (required)

        -o, --output
                    Output directory (required)

        -a, --aln-on-asm
                    Reads aligned on assembly (SAM format) (not recommended)

        -q, --output-read-groups
                    Output read groups (txt format)

        -p, --polish
                    Use this option if the input assembly is not polished

        -s, --dont_simplify
                    Don't rename the contigs and don't merge them

        --path-to-minimap2
                    Path to the executable minimap2 (if not in PATH)

        --path-to-racon
                    Path to the executable racon (if not in PATH)

        -F, --force Force overwrite of output folder if it exists
        -t, --threads
                    Number of threads
```

# Issues
 Most installation issues that we have seen yet stem from the use of too old compilers. Hairsplitter has been developed using gcc=11.2.0. Sometimes the default version of the compiler is too old (especially on servers). Specify gcc versions manually to cmake using `-DCMAKE_CXX_COMPILER=/path/to/modern/g++` and `-DCMAKE_C_COMPILER=/path/to/modern/gcc`.




