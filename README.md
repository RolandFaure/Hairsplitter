# Hairsplitter

Splits contigs into their different haplotypes (or repeats into their different versions).

# What is Hairsplitter ?

`Hairsplitter` takes as input an assembly (obtained by any means) and the long reads (including high-error rate long reads) used to build this assembly. For each contig it checks if the contig was built using reads from different haplotypes/regions. If it was, `Hairsplitter` separates the reads into as many groups as necessary and computes the different versions (e.g. alleles) of the contig actually present in the genome. It outputs a new assembly, where different versions of contigs are not collapsed into one but assembled separately.

# Why is it useful ?

`Hairsplitter` is mainly useful if you are trying to obtain a phased assembly. The main advantage of `Hairsplitter` compared to other techniques is that it is totally parameters-free. Most importantly, it does not requires to know the ploidy of the organism, and can infer different ploidies corresponding to different contigs. It can thus be used just as well on haploid assemblies (to improve the assembly of duplications) as on complex allotetraploids (to assemble separately the haplotypes). Just run the assembly through!

# Installation

A conda package is in preparation but is not available yet. For now, it is necessary to download and compile the code.

## Dependancies

- [minimap2](https://github.com/lh3/minimap2): this is not absolutely necessary. Another aligner can be used to manually align the reads.
- [racon](https://github.com/isovic/racon)
- CMake >= 3.8.12, make, C++11
- [GraphUnzip](https://github.com/nadegeguiglielmoni/GraphUnzip): this is included in the `Hairsplitter` folder in Hairsplitter/GraphUnzip/graphunzip.py, you do not need to install it separately. It is a python file, thus it will not work if taken out of its folder.
 
## Download & Compilation

To download and compile, run
```
git clone https://github.com/RolandFaure/Hairsplitter.git
cd Hairsplitter
cmake .
make
```

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
                    Reads aligned on assembly (PAF format)

        -q, --outputGAF
                    Output GAF file

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




