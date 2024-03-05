# GenomeTailor

## Installation

### Dependencies

To compile the project you will need CMake>=3.8.12 and GCC (tested with GCC=11.3.1).
You will need to have `racon`, `minimap2` and optionnally `minigraph`.

### Download & install

```
git clone https://github.com/RolandFaure/GenomeTailor.git
cd GenomeTailor
mkdir build && cd build
cmake ..
make
```

## Usage

```
SYNOPSIS
        build/GenomeTailor -i <input_assembly> -r <input_reads> -o <output_assembly> [-g <gaf_file>]
                           [-n <minigraph>] [-m <minimap2>] [-r <racon>] [-h] [-v]

OPTIONS
        -i, --input_assembly
                    input assembly in gfa format

        -r, --input_reads
                    input reads in fasta/q format

        -o, --output_assembly
                    output assembly in gfa format

        -g, --gaf_file
                    gaf file. Will be generated with minigraph if not provided

        -n, --minigraph
                    path to minigraph

        -m, --minimap2
                    path to minimap2

        -r, --racon path to racon
        -h, --help  print this help message and exit
        -v, --version
                    print version information and exit
```
