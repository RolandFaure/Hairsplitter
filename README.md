# Hairsplitter

Splits contigs into their different haplotypes (or repeats into their different versions).

# What is Hairsplitter ?

`Hairsplitter` takes as input an assembly (obtained by any means) and the long reads (including high-error rate long reads) used to build this assembly. For each contig it checks if the contig was built using reads from different haplotypes/regions. If it was, `Hairsplitter` separates the reads into as many groups as necessary and computes the different versions (e.g. alleles) of the contig actually present in the genome. It outputs a new assembly, where different versions of contigs are not collapsed into one but assembled separately.

# Why is it useful ?

`Hairsplitter` can be used to refine a metagenomic assembly. Assemblers commonly collapse closely related strains as on single genome. HairSplitter can recover the lost strains. The uncollapsed parts of the assembly are left as is.
`HairSplitter` is also useful for single-organism assembly, especially if you are trying to obtain a phased assembly. The main advantage of `Hairsplitter` compared to other techniques is that it is totally parameter-free. Most importantly, it does not requires to know the ploidy of the organism, and can infer different ploidies corresponding to different contigs. It can thus be used just as well on haploid assemblies (to improve the assembly of duplications) as on complex allotetraploids (to assemble separately the haplotypes). Just run the assembly through!

# Installation

A conda package is in preparation but is not available yet. For now, it is necessary to download and compile the code.

## Dependancies

### Quick conda dependancies

You can create and activate a conda environment with all dependencies installed by typing: 
```
conda create -c bioconda -c conda-forge -c anaconda -n hairsplitter minimap2 racon samtools cmake gxx gcc python scipy numpy
conda activate hairsplitter
```

### List of dependencies

- [minimap2](https://github.com/lh3/minimap2)
- [racon](https://github.com/isovic/racon)
- [samtools](www.htslib.org)
- CMake >= 3.8.12, make, gcc >= 11, g++ >= 11
- Python3 with numpy and scipy

If Minimap2, Racon or samtools are not in the PATH, their location should be specified through the `--path-to-minimap2`, `--path-to-racon` or `--path-to-samtools` options.
 
## Download & Compilation

To download and compile, run
```
git clone https://github.com/RolandFaure/Hairsplitter.git
cd Hairsplitter/src
mkdir build && cd build
cmake ..
make
```

# Usage

## Quick start

Let's say `reads.fastq` (ONT reads) were used to build assembly `assembly.gfa` (with any assembler)(the assembly can be in gfa or fasta format). To improve/phase the assembly using `Hairsplitter`, run
```
python hairsplitter.py -f reads.fastq -i assembly.gfa -x ont -o hairsplitter_out/
```

In the folder hairsplitter_out, you will find the new assembly, named `hairsplitter\_assembly.gfa`. Another generated file is `hairsplitter\_summary.txt`, in which are written which contigs are duplicated and merged.

## Options

```bash
usage: hairsplitter.py [-h] -i ASSEMBLY -f FASTQ [-x TECHNOLOGY] [-t THREADS] [-s] -o OUTPUT [-F] [--path_to_minimap2 PATH_TO_MINIMAP2]
                       [--path_to_racon PATH_TO_RACON] [--path_to_samtools PATH_TO_SAMTOOLS] [-d] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i ASSEMBLY, --assembly ASSEMBLY
                        Original assembly in GFA or FASTA format (required)
  -f FASTQ, --fastq FASTQ
                        Sequencing reads fasta/q (required)
  -x TECHNOLOGY, --technology TECHNOLOGY
                        {ont, pacbio, hifi} [ont]
  -t THREADS, --threads THREADS
                        Number of threads [1]
  -s, --dont_simplify   Don't rename the contigs and don't merge them
  -o OUTPUT, --output OUTPUT
                        Output directory
  -F, --force           Force overwrite of output folder if it exists
  --path_to_minimap2 PATH_TO_MINIMAP2
                        Path to the executable minimap2 [minimap2]
  --path_to_racon PATH_TO_RACON
                        Path to the executable racon [racon]
  --path_to_samtools PATH_TO_SAMTOOLS
                        Path to samtools [samtools]
  -d, --debug           Debug mode
  -v, --version         Print version and exit

```

# Issues
 Most installation issues that we have seen yet stem from the use of too old compilers. Hairsplitter has been developed using gcc=11.2.0. Sometimes the default version of the compiler is too old (especially on servers). Specify gcc versions manually to cmake using `-DCMAKE_CXX_COMPILER=/path/to/modern/g++` and `-DCMAKE_C_COMPILER=/path/to/modern/gcc`.




