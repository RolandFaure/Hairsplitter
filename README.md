![HS_logo](HS_logo.png)

Splits contigs into their different haplotypes (or repeats into their different versions).

*For developers working on similar problems:* HairSplitter is puposefully built as a series of modules that could be integrated in other software. See the [How does it work section](#work) and do not hesitate to contact me.

# What is Hairsplitter ?

`Hairsplitter` takes as input an assembly (obtained by any means) and the long reads (including high-error rate long reads) used to build this assembly. For each contig it checks if the contig was built using reads from different haplotypes/regions. If it was, `Hairsplitter` separates the reads into as many groups as necessary and computes the different versions (e.g. alleles) of the contig actually present in the genome. It outputs a new assembly, where different versions of contigs are not collapsed into one but assembled separately.

# Why is it useful ?

`Hairsplitter` can be used to refine a metagenomic assembly. Assemblers commonly collapse closely related strains as on single genome. HairSplitter can recover the lost strains. The uncollapsed parts of the assembly are left as is.
`HairSplitter` is also useful for single-organism assembly, especially if you are trying to obtain a phased assembly. The main advantage of `Hairsplitter` compared to other techniques is that it is totally parameter-free. Most importantly, it does not requires to know the ploidy of the organism, and can infer different ploidies corresponding to different contigs. It can thus be used just as well on haploid assemblies (to improve the assembly of duplications) as on complex allotetraploids (to assemble separately the haplotypes). Just run the assembly through!

# Installation

A conda package is in preparation but is not available yet. For now, it is necessary to download and compile the code.

## Dependancies

### Quick conda dependancies

The recommended way to install HairSplitter is to create and activate a conda/mamba environment with all dependencies: 
```
conda create -c bioconda -c conda-forge -c anaconda -n hairsplitter cmake gxx gcc python scipy numpy minimap2 racon htslib==1.16 samtools medaka
conda activate hairsplitter
```

### List of dependencies

- [minimap2](https://github.com/lh3/minimap2)
- [racon](https://github.com/isovic/racon) and/or [medaka](https://github.com/nanoporetech/medaka)
- [samtools](www.htslib.org)
- CMake >= 3.8.12, make, gcc >= 11, g++ >= 11
- Python3 with numpy and scipy

If Minimap2, Racon, Medaka or samtools are not in the PATH, their location should be specified through the `--path-to-minimap2`, `--path-to-racon`, `path-to-medaka` or `--path-to-samtools` options.
 
## Download & Compilation

To download and compile, run
```
git clone https://github.com/RolandFaure/Hairsplitter.git
cd Hairsplitter/src/
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
  -m, --multiploid      Use this option if all haplotypes can be assumed to have the same coverage
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
 
 <a name="work">
</a>

# How does it work ?

HairSplitter is organized as series of modules, some of these modules being of independant interest. The full documentation can be found in the doc/ folder.

1. *Cleaning the assembly* Ideally, the assembly would be purged of all assembly errors. In practice, ensure there is no over-duplication by deleting unconnected contigs that align very well on other contigs.

2. *Calling variants* Variants are called using an alignment of the reads on the assembly. For now, a basic pileup is used. Calling variants in a metagenomic context is hard: favor calling false variants over missing true variants - the false variants will be filtered afterward.

3. *Filtering variants* This step is crucial. Each called variant partition the reads in groups. Keep only variants which partition occur frequently, because this cannot be chance. This way, only very robust variant are kept.

4. *Separating the reads* Based on the robust variants, HairSplitter inspect each contig and determine if several distinct groups of reads align there. If it is the case, it means that several different versions of the contig exist.

5. *Creating the new contigs* Create every new contig by polishing the existing contig using the several groups of reads.

6. *Improving contiguity* Contigs are generally separated only locally. To improve contiguity, use the long reads that align on several contigs sequentially.
 
# Citation
 A preprint will very soon be available online. Until then, you can cite the [proceedings of the JOBIM conference](jobim2023.sciencesconf.org/data/pages/proceedings.pdf), page 124:
```
HairSplitter: separating strains in metagenome assemblies with long reads
Roland Faure, Jean-Francois Flot, Dominique Lavenier
JOBIM 2023 proceedings, p.124-131
```
 
 




