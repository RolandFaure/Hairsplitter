# HairSplitter

This files contains the description of the different modules of HairSplitter and the specific file formats used for each module. To use HairSplitter, refer to the README at the root of the folder (`github.com/RolandFaure/HairSplitter`).
Do not hesitate to use the issues of the github to suggest improvements or new features in the way the modules are designed.

## Overview

HairSplitter is organized as series of modules, some of these modules being of independant interest. 

1. *Cleaning the assembly* Ideally, the assembly would be purged of all assembly errors. In practice, ensure there is no over-duplication by deleting unconnected contigs that align very well on other contigs.

2. *Calling variants* Variants are called using an alignment of the reads on the assembly. For now, a basic pileup is used. Calling variants in a metagenomic context is hard: favor calling false variants over missing true variants - the false variants will be filtered afterward.

3. *Filtering variants* This step is crucial. Each called variant partition the reads in groups. Keep only variants which partition occur frequently, because this cannot be chance. This way, only very robust variant are kept.

4. *Separating the reads* Based on the robust variants, HairSplitter inspect each contig and determine if several distinct groups of reads align there. If it is the case, it means that several different versions of the contig exist.

5. *Creating the new contigs* Create every new contig by polishing the existing contig using the several groups of reads.

6. *Improving contiguity* Contigs are generally separated only locally. To improve contiguity, use the long reads that align on several contigs sequentially.



