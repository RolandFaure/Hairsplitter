# HairSplitter

This files contains the description of the different modules of HairSplitter and the specific file formats used for each module. To use HairSplitter, refer to the README at the root of the folder (`github.com/RolandFaure/HairSplitter`).
Do not hesitate to use the issues of the github to suggest improvements or new features in the way the modules and/or the formats are designed.

## Overview

HairSplitter is organized as series of modules, some of these modules being of independant interest. 

1. *Cleaning the assembly*. Ideally, the assembly would be purged of all assembly errors. In practice, ensure there is no over-duplication by deleting unconnected contigs that align very well on other contigs.

2. *Calling variants*. Variants are called using an alignment of the reads on the assembly. For now, a basic pileup is used. Calling variants in a metagenomic context is hard: favor calling false variants over missing true variants - the false variants will be filtered afterward.

3. *Filtering variants*. This step is crucial. Each called variant partition the reads in groups. Keep only variants which partition occur frequently, because this cannot be chance. This way, only very robust variant are kept.

4. *Separating the reads*. Based on the robust variants, HairSplitter inspect each contig and determine if several distinct groups of reads align there. If it is the case, it means that several different versions of the contig exist.

5. *Creating the new contigs*. Create every new contig by polishing the existing contig using the several groups of reads.

6. *Improving contiguity*. Contigs are generally separated only locally. To improve contiguity, use the long reads that align on several contigs sequentially.

## Module: cleaning the assembly

```
Usage: clean_graph <assemblyFile> <newAssembly> <outputFolder> <logFile> <num_threads> <path_to_minimap>
```
**assemblyFile** Original assembly in GFA format

**newAssembly** Output, in GFA format

**outputFolder** Folder to store temporary file (typically alignment files). Can be e.g. "hairsplitter/tmp".

**logFile** File to write what clean_graph does

**num_threads** Number of threads

**path_to_minimap** Path to executable minimap2

## Module: Calling variants

```
Usage: ./call_variants <gfa_file> <reads_file> <sam_file> <num_threads> <tmpDir> <error_rate_out> <DEBUG> <columns_out>
```

**gfa_file** Assembly in GFA format

**reads_file** Reads in fasta (.fa or .fasta) or fastq (.fq or .fastq) format

**sam_file** Alignment of the reads on the assembly. SAM format.

**num_threads** Number of threads

**tmpDir** Folder to store temporary file. Can be e.g. "hairsplitter/tmp".

**error_rate_out** File where the computed error rate of the reads will be stored. e.g. "hairsplitter/tmp/error_rate.txt". `call\_variants` will fill it with a single floating point value.

**DEBUG** 0 or 1. 

**columns_out** File containing the position of variants and the value of each read at each position in [COL format](#col).

## Module: Filtering variants

```
Usage: ./filter_variants <columns> <error_rate> <num_threads> <DEBUG> <columns_out>
```

**columns** File containing the position of variants and the value of each read at each position in [COL format](#col).

**error_rate** File containing the error rate of the reads. Plain file containing single floating point value.

**num_threads** Number of threads

**DEBUG** 0 or 1.

**columns_out** File containing the position of robust variants and the value of each read at each position in [COL format](#col).

## Module: Separating the reads

```
./separate_reads <columns> <num_threads> <error_rate> <DEBUG> <outfile>
```

**columns** File containing the position of robust variants and the value of each read at each position in [COL format](#col).

**num_threads** Number of threads

**error_rate** File containing the error rate of the reads. Plain file containing single floating point value.

**DEBUG** 0 or 1.

**outfile** File containing the positions on the assembly where several groups of reads are present, and the composition of the groups. In [GRO format](#gro)

## Module: Creating new contigs

```
Usage: ./create_new_contigs <assembly> <reads_file> <error_rate> <gro_file> <tmpfolder> <num_threads> <technology> <output_graph> <output_gaf> <path_to_minimap> <path-to-racon> <debug>
```
**assembly** Cleaned assembly in GFA format

**reads_file** Reads in fasta (.fa or .fasta) or fastq (.fq or .fastq) format

**error_rate** File containing the error rate of the reads. Plain file containing single floating point value.

**gro_file** File containing the positions on the assembly where several groups of reads are present, and the composition of the groups. In [GRO format](#gro)

**tmpfolder** Folder to store temporary file. Can be e.g. "hairsplitter/tmp".

**num_threads** Number of threads

**output_gaf** File created by create_new_contigs that describes how the reads align on the contig successively in the [GAF format](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf). Only the first 6 fields need to be filled 

**path_to_minimap** Path to executable minimap2

**path_to_racon** Path to executable racon

**debug** 0 or 1.

## Module: Improving contiguity

HairSplitter uses [GraphUnzip](github.com/nadegeguiglielmoni/GraphUnzip) for this part. Please refer to this page.

<a name="col">
</a>

## File format specification: COL
The goal of the format is to inventoriate the position of variants and which reads contain which variants.
It is build as a tab-separated file with the first column being the type of the line. There are 3 types of lines:
### Contig lines
A contig line describes a contig of the assembly.
```
CONTIG  edge_1  6793    6.84705
```
4 fields in this kind of line:
1. CONTIG
2. Name of contig (/!\ no tab in the name)
3. Length of contig
4. Coverage depth of the contig (float)

### Read lines
Read lines follow a contig line and describe all the reasd that aling on the contig.
```
READ    @0_7b15b848-cb53-91a7-b204-e2cced890836      0       3660    148     3809    1
```
7 fields in this kind of line:
1. READ
2. Name of read (/!\ no tab in the name)
3. Coordinates of the beginning of the alignment on the read
4. Coordinates of the end of the alignment on the read
5. Coordinates of the beginning of the alignment on the contig
6. Coordinates of the end of the alignment on the contig
7. Strand of alignment (0: reverse, 1: forward)

### SNPS lines
SNPS line come after the read lines and describe the position of the pileup of reads on the contig.
```
SNPS    1259       A       T       :AAAA TCAT TTTCAA TTT  TTTTAAAAAAA AATATATATTTAAAA 
```
5 fields in this kind of line:
1. SNPS
2. Position of the SNP on the contig
3. Majority allele at this position (can be any ASCII character)
4. Minority allele at this position (can be any ASCII character)
5. Pileup. Must starts by ":". Then, there must be exactly as many character as there are reads aligning on a contig. The character at position _n_ describes the allele of the _n_th read (corresponding to the _n_th READ line after the last CONTIG line) at this position. ' ' is inserted if the read is not defined at this position.

A complete (small) COL file can look like this:
```
CONTIG	edge_1	1034	8.3
READ	read_1	3300	4356	0	1034	1
READ	read_4	0	1056	0	1034	0
READ	read_5	0	567	345	901	1
READ	read_6	100	1203	0	1034	1
SNPS	23	A	T	:  AT
SNPS	345	C	-	:CC--
SNPS	678	-	G	:-G -
```

<a name="GRO">
</a>

## File format specification: GRO
The goal of the format is to inventoriate the groups of reads along the contigs of the assembly.
It is build as a tab-separated file with the first column being the type of the line. There are 3 types of lines:
### Contig lines
A contig line describes a contig of the assembly.
```
CONTIG  edge_1  6793    6.84705
```
4 fields in this kind of line:
1. CONTIG
2. Name of contig (/!\ no tab in the name)
3. Length of contig
4. Coverage depth of the contig (float)

### Read lines
Read lines follow a contig line and describe all the reasd that aling on the contig.
```
READ    @0_7b15b848-cb53-91a7-b204-e2cced890836      0       3660    148     3809    1
```
7 fields in this kind of line:
1. READ
2. Name of read (/!\ no tab in the name)
3. Coordinates of the beginning of the alignment on the read
4. Coordinates of the end of the alignment on the read
5. Coordinates of the beginning of the alignment on the contig
6. Coordinates of the end of the alignment on the contig
7. Strand of alignment (0: reverse, 1: forward)

### Group lines
Group lines follow read lines and describe how reads are split along the contig.
```
GROUP   6000    7999    -2,0,0,-2,0,-2,0,-2,-2,-2,0,-2,-2,1,-2,1,-2,1,1,-2,-2,
```
4 fields in this kind of line:
1. GROUP
2. Start position of the group on the contig
3. End position of the group on the contig
4. Assignment of the reads in groups. Comma-separated list of integers. The integer at position _n_ is the identifier of the group of the _n_th read (corresponding to the _n_th READ line after the last CONTIG line). `-2` is a special int to describe that the read is not defined on this interval. In the above example, there is thus only two groups of reads.

The GROUP lines MUST tile the entire contig!

An example of a GRO file is:
```
CONTIG	edge_1	1034	8.3
READ	read_1	3300	4356	0	1034	1
READ	read_4	0	1056	0	1034	0
READ	read_5	0	567	345	901	1
READ	read_6	100	1203	0	1034	1
GROUP	0	199	0,0,-2,0
GROUP	200	344	0,0,1,1
GROUP	345	901	4,4,4,4
GROUP	902	1034	0,0,-2,2
GROUP	1034	1034	-2,-2,-2,-2
```





