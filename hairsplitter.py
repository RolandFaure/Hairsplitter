#!/usr/bin/env python3

"""
HairSplitter takes as input a collapsed assembly and a set of long reads and outputs an uncollapsed assembly.
This is the master file that calls all the other scripts.
Author: Roland Faure

"""

__author__ = "Roland Faure"
__license__ = "GPL3"
__version__ = "1.8.2"
__date__ = "2024-04-30"
__maintainer__ = "Roland Faure"
__email__ = "roland.faure@irisa.fr"
__github__ = "github.com/RolandFaure/HairSplitter"
__status__ = "Prototype"

import sys
import os
import argparse
import datetime


# command line arguments
def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--assembly", help="Original assembly in GFA or FASTA format (required)", required=True)
    parser.add_argument("-f", "--fastq", help="Sequencing reads fastq or fasta (required)", required=True)
    parser.add_argument("-c", "--haploid-coverage", help="Expected haploid coverage. 0 if does not apply [0]", default=0)
    parser.add_argument("-x", "--technology", help="{ont, pacbio, hifi} [ont]", default="ont")
    parser.add_argument("-p", "--polisher", help="{racon,medaka} medaka is more accurate but much slower [racon]", default="racon")
    # parser.add_argument("-m", "--multiploid", help="Use this option if all haplotypes can be assumed to have the same coverage", action="store_true")
    parser.add_argument("--correct-assembly", help="Correct structural errors in the input assembly (time-consuming)", action="store_true")
    parser.add_argument("-t", "--threads", help="Number of threads [1]", default=1)
    parser.add_argument("-o", "--output", help="Output directory", required=True)
    parser.add_argument("--resume", help="Resume from a previous run", action="store_true")
    parser.add_argument("-s", "--dont_simplify", help="Don't merge the contig", action="store_true")
    parser.add_argument("-P", "--polish-everything", help="Polish every contig with racon, even those where there is only one haplotype ", action="store_true")
    parser.add_argument("-F", "--force", help="Force overwrite of output folder if it exists", action="store_true")
    parser.add_argument("-l", "--low-memory", help= "Turn on the low-memory mode (at the expense of speed)", action="store_true")
    parser.add_argument("--clean", help="Clean the temporary files", action="store_true")
    parser.add_argument("--rarest-strain-abundance", help="Limit on the relative abundance of the rarest strain to detect (0 might be slow for some datasets) [0.01]", default=0.01, type=float)
    parser.add_argument("--minimap2-params", help="Parameters to pass to minimap2", default="")
    parser.add_argument("--path_to_minimap2", help="Path to the executable minimap2 [minimap2]", default="minimap2")
    parser.add_argument("--path_to_minigraph", help="Path to the executable minigraph [minigraph]", default="minigraph")
    parser.add_argument("--path_to_racon", help="Path to the executable racon [racon]", default="racon")
    parser.add_argument("--path_to_medaka", help="Path to the executable medaka [medaka]", default="medaka")
    parser.add_argument("--path_to_samtools", help="Path to samtools [samtools]", default="samtools")
    parser.add_argument("--path_to_python", help="Path to python [python]", default="python")
    parser.add_argument("--path_to_raven", help="Path to raven [raven]", default="raven")

    parser.add_argument("-v", "--version", help="Print version and exit", action="store_true")
    parser.add_argument("-d", "--debug", help="Debug mode", action="store_true")


    return parser.parse_args()

def check_dependencies(tmp_dir, minimap2, minigraph, racon, medaka, polisher, samtools, path_to_src, path_to_python, skip_minigraph\
                       , path_GenomeTailor, path_cut_gfa, path_fa2gfa, path_gfa2fa, path_call_variants, path_separate_reads, path_create_new_contigs\
                        , path_determine_multiplicity,  path_graphunzip\
                        , path_raven):

    com = " --version > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
    mini_run = os.system(minimap2 + com)

    minigraph_run = os.system(minigraph + com)

    if polisher != "medaka" :
        racon_run = os.system(racon + com)

    if polisher != "racon" :
        medaka_run = os.system(medaka + com)

    if not skip_minigraph :
        raven_run = os.system(path_raven + com)

    samtools_run = os.system(samtools + com)

    command = path_to_python + " --version > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
    python_run = os.system(command)

    #print a table listing the dependencies that are ok or not
    print("\n===== Checking dependencies =====\n")
    print("______________________________________________________________")
    print("|  Dependency  |  Status  |            Path Tried            |")
    print("|--------------|----------|----------------------------------|")
    #if the dependancy is ok, print OK in green, else print ERROR in red
    if mini_run == 0:
        print("| minimap2     |   \033[92mOK\033[0m     | "+minimap2, end="")
        #add white spaces to align the columns
        for i in range(0, 33-len(minimap2)):
            print(" ", end="")
        print("|")
    else:
        print("| minimap2     |  \033[91mERROR\033[0m   | "+minimap2, end="")
        for i in range(0, 33-len(minimap2)):
            print(" ", end="")
        print("|")

    if minigraph_run == 0 and not skip_minigraph:
        print("| minigraph    |   \033[92mOK\033[0m     | "+minigraph, end="")
        for i in range(0, 33-len(minigraph)):
            print(" ", end="")
        print("|")
    elif not skip_minigraph:
        print("| minigraph    |  \033[91mERROR\033[0m   | "+minigraph, end="")
        for i in range(0, 33-len(minigraph)):
            print(" ", end="")
        print("|")

    if polisher != "medaka" and racon_run == 0:
        print("| racon        |   \033[92mOK\033[0m     | "+racon, end="")
        for i in range(0, 33-len(racon)):
            print(" ", end="")
        print("|")
    elif polisher != "medaka":
        print("| racon        |  \033[91mERROR\033[0m   | "+racon, end="")
        for i in range(0, 33-len(racon)):
            print(" ", end="")
        print("|")

    if polisher != "racon" and medaka_run == 0:
        print("| medaka       |   \033[92mOK\033[0m     | "+medaka, end="")
        for i in range(0, 33-len(medaka)):
            print(" ", end="")
        print("|")
    elif polisher != "racon":
        print("| medaka       |  \033[91mERROR\033[0m   | "+medaka, end="")
        for i in range(0, 33-len(medaka)):
            print(" ", end="")
        print("|")

    if not skip_minigraph and raven_run == 0:
        print("| raven        |   \033[92mOK\033[0m     | "+path_raven, end="")
        for i in range(0, 33-len(path_raven)):
            print(" ", end="")
        print("|")
    elif not skip_minigraph:
        print("| raven        |  \033[91mERROR\033[0m   | "+path_raven, end="")
        for i in range(0, 33-len(path_raven)):
            print(" ", end="")
        print("|")

    if samtools_run == 0:
        print("| samtools     |   \033[92mOK\033[0m     | "+samtools, end="")
        for i in range(0, 33-len(samtools)):
            print(" ", end="")
        print("|")
    else:
        print("| samtools     |  \033[91mERROR\033[0m   | "+samtools, end="")
        for i in range(0, 33-len(samtools)):
            print(" ", end="")
        print("|")

    if python_run == 0:
        print("| python       |   \033[92mOK\033[0m     | "+path_to_python, end="")
        for i in range(0, 33-len(path_to_python)):
            print(" ", end="")
        print("|")
    else:
        print("| python       |  \033[91mERROR\033[0m   | "+path_to_python, end="")
        for i in range(0, 33-len(path_to_python)):
            print(" ", end="")
        print("|")
    
    print("______________________________________________________________\n")

    #if any of the dependencies is not ok, exit
    if mini_run != 0 or (minigraph_run != 0 and not skip_minigraph) or (polisher != "medaka" and racon_run != 0) or (polisher != "racon" and medaka_run != 0) or samtools_run != 0 or python_run != 0 \
        or (not skip_minigraph and raven_run != 0):
        print("ERROR: Some dependencies could not run. Check the path to the executables.")
        sys.exit(1)

    #now check the path to our own executable, to make sure the installation ran smoothly
    command = path_GenomeTailor + " --version > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
    GenomeTailor_run = os.system(command)
    if GenomeTailor_run != 0:
        command = "HS_GenomeTailor --version > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
        GenomeTailor_run = os.system(command)
        if GenomeTailor_run != 0:
            print("ERROR: GenomeTailor could not run. Problem in the installation.")
            print("Was trying to run first: " + path_GenomeTailor + " --version > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt")
            print("Was trying to run: " + command)
            sys.exit(1)
        else:
            path_GenomeTailor = "HS_GenomeTailor"

    command = path_cut_gfa + " -h > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
    cut_gfa_run = os.system(command)
    if cut_gfa_run != 0:
        command = "cut_gfa.py -h > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
        cut_gfa_run = os.system(command)
        if cut_gfa_run != 0:
            print("ERROR: cut_gfa.py could not run. Problem in the installation.")
            print("Was trying to run first: " + path_cut_gfa + " --help > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt")
            print("Was trying to run: " + command)
            sys.exit(1)
        else:  
            path_cut_gfa = "cut_gfa.py"

    command = path_fa2gfa + " --version > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
    fa2gfa_run = os.system(command)
    if fa2gfa_run != 0:
        command = "HS_fa2gfa --version > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
        fa2gfa_run = os.system(command)
        if fa2gfa_run != 0:
            print("ERROR: fa2gfa could not run. Problem in the installation.")
            print("Was trying to run: " + command)
            sys.exit(1)
        else:
            path_fa2gfa = "HS_fa2gfa"

    command = path_gfa2fa + " --version > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
    gfa2fa_run = os.system(command)
    if gfa2fa_run != 0:
        command = "HS_gfa2fa --version > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
        gfa2fa_run = os.system(command)
        if gfa2fa_run != 0:
            print("ERROR: gfa2fa could not run. Problem in the installation.")
            print("Was trying to run: " + command)
            sys.exit(1)
        else:
            path_gfa2fa = "HS_gfa2fa"

    command = path_call_variants + " --version > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
    call_variants_run = os.system(command)
    if call_variants_run != 0:
        command = "HS_call_variants --version > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
        call_variants_run = os.system(command)
        if call_variants_run != 0:
            print("ERROR: call_variants could not run. Problem in the installation.")
            print("Was trying to run: " + command)
            sys.exit(1)
        else:
            path_call_variants = "HS_call_variants"

    command = path_separate_reads + " --help > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
    separate_reads_run = os.system(command)
    if separate_reads_run != 0:
        command = "HS_separate_reads --help > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
        separate_reads_run = os.system(command)
        if separate_reads_run != 0:
            print("ERROR: separate_reads could not run. Problem in the installation.")
            print("Was trying to run first: " + path_separate_reads + " --help > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt")
            print("Was trying to run: " + command)
            sys.exit(1)
        else:
            path_separate_reads = "HS_separate_reads"

    command = path_create_new_contigs + " --help > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
    create_new_contigs_run = os.system(command)
    if create_new_contigs_run != 0:
        command = "HS_create_new_contigs --help > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
        create_new_contigs_run = os.system(command)
        if create_new_contigs_run != 0:
            print("ERROR: create_new_contigs could not run. Problem in the installation.")
            print("Was trying to run: " + command)
            sys.exit(1)
        else:
            path_create_new_contigs = "HS_create_new_contigs"

    command = path_graphunzip + " unzip --help > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
    graphunzip_run = os.system(command)
    if graphunzip_run != 0:
        command = "graphunzip.py unzip --help > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
        graphunzip_run = os.system(command)
        if graphunzip_run != 0:
            print("ERROR: graphunzip.py could not run. Problem in the installation.")
            print("Was trying to run first: " + path_graphunzip + " unzip --help > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt")
            print("Was trying to run: " + command)
            sys.exit(1)
        else:
            path_graphunzip = "graphunzip.py"

    command = path_determine_multiplicity + " --help > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
    determine_multiplicity_run = os.system(command)
    if determine_multiplicity_run != 0:
        command = "determine_multiplicity.py --help > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
        determine_multiplicity_run = os.system(command)
        if determine_multiplicity_run != 0:
            print("ERROR: determine_multiplicity.py could not run. Problem in the installation.")
            print("Was trying to run first: " + path_determine_multiplicity + " --help > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt")
            print("Was trying to run: " + command)
            sys.exit(1)
        else:
            path_determine_multiplicity = "determine_multiplicity.py"

    return path_GenomeTailor, path_cut_gfa, path_fa2gfa, path_gfa2fa, path_call_variants, path_separate_reads, path_create_new_contigs, path_determine_multiplicity, path_graphunzip

def main():

    print("\n\t******************\n\t*                *\n\t*  Hairsplitter  *\n\t*    Welcome!    *\n\t*                *\n\t******************\n")
    sys.stdout.flush()

    if len(sys.argv) > 1 and (sys.argv[1] == "-v" or sys.argv[1] == "--version"):
        print("HairSplitter v"+__version__+" ("+__github__+"). Last update: "+__date__)
        sys.exit(0)

    args = parse_args()
    nb_threads = args.threads
    #path to src folder can be imputed from the first argument of the command line
    path_to_src = sys.argv[0].split("hairsplitter.py")[0]+"src/"
    path_to_minimap2 = args.path_to_minimap2
    path_to_minigraph = args.path_to_minigraph
    path_to_raven = args.path_to_raven
    readsFile = args.fastq
    tmp_dir = args.output.rstrip('/') + "/tmp"
    path_to_python = args.path_to_python
    low_memory = args.low_memory
    skip_minigraph = not args.correct_assembly
    rarest_strain_abundance = args.rarest_strain_abundance
    minimap2_params = args.minimap2_params
    haploid_coverage = float(args.haploid_coverage)
    continue_from_previous_run = args.resume
    clean_tmp = args.clean

    path_GenomeTailor = path_to_src + "build/HS_GenomeTailor/HS_GenomeTailor"
    path_cut_gfa = path_to_python + " " + path_to_src + "cut_gfa.py"
    path_fa2gfa = path_to_src + "build/HS_fa2gfa"
    path_gfa2fa = path_to_src + "build/HS_gfa2fa"
    path_call_variants = path_to_src + "build/HS_call_variants"
    path_separate_reads = path_to_src + "build/HS_separate_reads"
    path_create_new_contigs = path_to_src + "build/HS_create_new_contigs"
    path_graphunzip = path_to_python + " " +path_to_src + "GraphUnzip/graphunzip.py"
    path_determine_multiplicity = path_to_python + " " + path_to_src + "GraphUnzip/determine_multiplicity.py"

    logFile = args.output.rstrip('/') + "/hairsplitter.log"

    # check if output folder exists
    if os.path.exists(args.output) and not args.force and not args.resume:
        print("ERROR: output folder already exists. Use -F to overwrite.")
        sys.exit(1)
    elif not os.path.exists(args.output) :
        # create output folder
        os.mkdir(args.output)

    #output the command line used to run HairSplitter and the version in the log file
    f = open(logFile, "w")
    f.write("HairSplitter v"+__version__+" ("+__github__+"). Last update: "+__date__+"\n")
    f.write("Command line used to run HairSplitter: "+" ".join(sys.argv)+"\n")
    f.close()

    #print the command line used to run HairSplitter
    print(" ".join(sys.argv))
    print("HairSplitter v"+__version__+" ("+__github__+"). Last update: "+__date__)
    if args.version:
        sys.exit(0)

    polisher = args.polisher.lower()
    if polisher != "racon" and polisher != "medaka":
        print("ERROR: polisher must be either racon or medaka")
        f = open(logFile, "w")
        f.write("ERROR: polisher must be either racon or medaka\n")
        f.close()
        sys.exit(1)

    reads_on_asm = tmp_dir + "/reads_on_asm.sam"

    #check if all the files and dependencies are here

    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    # check if input files exist
    if not os.path.exists(args.assembly):
        print("ERROR: not found assembly (" + args.assembly + ")")
        sys.exit(1)
    if not os.path.exists(args.fastq):
        print("ERROR: not found fastq file (" + args.fastq + ")")
        sys.exit(1)

    elif not args.fastq.endswith(".fastq") and not args.fastq.endswith(".fq") and not args.fastq.endswith(".fastq.gz") and not args.fastq.endswith(".fq.gz") and not args.fastq.endswith(".fasta") and not args.fastq.endswith(".fa") and not args.fastq.endswith(".fna") and not args.fastq.endswith(".fasta.gz") and not args.fastq.endswith(".fa.gz"):
        print("ERROR: fastq file must be in FASTQ or FASTA format (potentially gzipped). File extension not recognized.")
        sys.exit(1)

    #check the dependencies
    path_GenomeTailor, path_cut_gfa, path_fa2gfa, path_gfa2fa, path_call_variants, path_separate_reads, path_create_new_contigs, path_determine_multiplicity, path_graphunzip =\
        check_dependencies(tmp_dir, args.path_to_minimap2, args.path_to_minigraph, args.path_to_racon, \
                       args.path_to_medaka, args.polisher, args.path_to_samtools, path_to_src, \
                        path_to_python, skip_minigraph, path_GenomeTailor, path_cut_gfa, path_fa2gfa, \
                        path_gfa2fa, path_call_variants, path_separate_reads, path_create_new_contigs, \
                        path_determine_multiplicity,
                        path_graphunzip, path_to_raven)

    #check the read file and unzip it if needed (converting it to fasta if in fastq)
    if readsFile[-3:] == ".gz":
        print("\n===== STAGE 0: Decompressing input reads [", datetime.datetime.now() ,"]\n\n")
        if continue_from_previous_run and os.path.exists(tmp_dir + "/reads.fasta") :
            print(" - Already decompressed reads file found from previous run")
            readsFile = tmp_dir + "/reads.fasta"
        else:
            continue_from_previous_run = False
            if readsFile[-6:-3] == ".fa" or readsFile[-9:-3] == ".fasta" :
                command = "gzip -d " + readsFile + " -c > " + tmp_dir + "/reads.fasta"
                readsFile = tmp_dir + "/reads.fasta"
            else :
                command = "gzip -d " + readsFile + " -c | sed -n '1~4s/^@/>/p;2~4p' > " + tmp_dir + "/reads.fasta"
                readsFile = tmp_dir + "/reads.fasta"
            print(" Running: " + command)
            res_gunzip = os.system(command)
            if res_gunzip != 0:
                print("ERROR: gzip failed. Was trying to run: " + command)
                sys.exit(1)
        
    # run the pipeline

    # 0. Convert the assembly to gfa if needed
    if args.assembly[-3:] == "gfa":
        gfaAssembly = args.assembly
    elif args.assembly[-5:] == "fasta" or args.assembly[-2:] == "fa" or args.assembly[-3:]=="fna":
        if not continue_from_previous_run or not os.path.exists(tmp_dir + "/assembly.gfa") :
            gfaAssembly = tmp_dir + "/assembly.gfa"
            command = path_fa2gfa + " " + args.assembly + " > " + gfaAssembly
            res_fasta2gfa = os.system(command)
            if res_fasta2gfa != 0:
                print("ERROR: Conversion from fasta to gfa failed while running the command:\n" + command)
                sys.exit(1)
    else:
        print("ERROR: Assembly file must be in GFA or FASTA format. File extension not recognized.")
        sys.exit(1)

    # 1. Clean the assembly using correct_structural_errors.py
    print("\n===== STAGE 1: Cleaning graph of hidden structural variations [", datetime.datetime.now() ,"]\n\n")
    print(" When several haplotypes are present, it is common that big structural variations between the haplotypes go unnoticed by the assembler. Here, HairSplitter corrects the assembly"
        " by making sure that all reads align end-to-end of the assembly. This module is now a standalone tool available at github.com/rolandfaure/genometailor\n")
    sys.stdout.flush()
    
    new_assembly = tmp_dir + "/cleaned_assembly.gfa"
    N50 = 0
    if not skip_minigraph :

        if continue_from_previous_run and os.path.exists(new_assembly) :
            print(" - Already cleaned assembly found from previous run")
        else:
            continue_from_previous_run = False
            command = path_GenomeTailor + " -i " + gfaAssembly + " -o " + new_assembly + " -r " + readsFile + " -t " + str(nb_threads) \
                + " --minimap2 " + args.path_to_minimap2 + " --minigraph " + args.path_to_minigraph + " --racon " + args.path_to_racon + " --path-to-raven " + path_to_raven \
                + " > " + tmp_dir + "/logGenomeTailor.txt"
            # command = "python " + path_to_src + "GraphUnzip/correct_structural_errors.py -a " + gfaAssembly + " -o " + new_assembly + " -r " + readsFile + " -t " \
            #     + str(nb_threads) + " --minimap2 " + args.path_to_minimap2 + " --minigraph " + args.path_to_minigraph + " --racon " + args.path_to_racon \
            #     + " --folder " + tmp_dir
            print(" Running: ", command)
            #write in the log file where to look in case of error
            f = open(logFile, "w")
            f.write("==== STAGE 1: Cleaning graph of hidden structural variations   ["+str(datetime.datetime.now())+"]\n")
            f.write(command)
            f.close()
            res_clean = os.system(command)
            if res_clean != 0:
                print("ERROR: Cleaning the assembly failed. Was trying to run: " + command)
                sys.exit(1)

            print(" - Improved alignment of reads on assembly. The improved assembly is stored in " + new_assembly)

            #now check if the improved assembly is not too complicated, else fall back on the original assembly
            #the metric is : did the N50 fall below 10kb ?
            f = open(new_assembly, "r")
            contigs_length = []
            for line in f :
                if "S" == line[0] :
                    contigs_length.append(len(line.split("\t")[2]))
            f.close()
            contigs_length.sort(reverse=True)
            cumul = 0
            total_length = sum(contigs_length)
            for l in contigs_length :
                cumul += l
                if cumul > total_length/2 :
                    N50 = l
                    break

    if N50 < 10000 and not skip_minigraph:
        print(" - The improved assembly is too complicated, falling back on the original assembly")
        new_assembly = gfaAssembly
    elif skip_minigraph :
        f = open(logFile, "w")
        f.write("==== STAGE 1: Cleaning graph of hidden structural variations   ["+str(datetime.datetime.now())+"]\n")
        f.write(" - Skipping the assembly correction step because --correct-assembly was not used")
        f.close()
        print(" - Skipping the assembly correction step because --correct-assembly was not used")
        new_assembly = gfaAssembly

    # 2. Map the reads on the assembly
    print("\n===== STAGE 2: Aligning reads on the reference   [", datetime.datetime.now() ,"]\n")
    sys.stdout.flush()

    # 2.1. Cut the contigs in chunks of 300000bp to avoid memory issues
    print(" - Cutting the contigs in chunks of 300000bp to avoid memory issues")
    command = path_cut_gfa + " -a " + new_assembly + " -l 300000 -o " + tmp_dir + "/cut_assembly.gfa"
    #write in the log file the time at which the alignment starts
    f = open(logFile, "a")
    f.write("\n==== STAGE 2: Aligning reads on the reference   ["+str(datetime.datetime.now())+"]\n")
    f.write(" - Cutting the contigs in chunks of 300000bp to avoid memory issues\n")
    f.write(command)

    res_cut_gfa = os.system(command)
    if res_cut_gfa != 0 :
        print("ERROR: cut_gfa.py failed. Was trying to run: " + command)
        sys.exit(1)
    new_assembly = tmp_dir + "/cut_assembly.gfa"

    # 2.2 Convert the assembly in fasta format
    print(" - Converting the assembly in fasta format")
    fastaAsm = tmp_dir + "/cleaned_assembly.fasta"
    command = path_gfa2fa + " " + new_assembly + " > " + fastaAsm
    f.write(" - Converting the assembly in fasta format\n")
    f.write(command)
    f.close()
    res_gfa2fasta = os.system(command)
    if res_gfa2fasta != 0 :
        print("ERROR: gfa2fa failed UUE. Was trying to run: " + command)
        sys.exit(1)

    techno_flag = ""
    technology = args.technology.lower()
    if technology == "pacbio" or technology == "pb":
        techno_flag = "-x map-pb"
    elif technology == "hifi" :
        techno_flag = "-x map-hifi"
    else :
        techno_flag = "-x map-ont"
    
    # 2.3 Align the reads on the assembly
    if not continue_from_previous_run or not os.path.exists(reads_on_asm) :
        continue_from_previous_run = False
        print(" - Aligning the reads on the assembly")

        #run minimap but do not store the sequences, they are still in the file of reads
        command = path_to_minimap2 + " " + fastaAsm + " " + readsFile + " " + techno_flag + " -a --secondary=no -M 0.05 -Y -t "+ str(nb_threads) + " " + minimap2_params + " 2> "+tmp_dir+"/logminimap.txt > " + reads_on_asm + " 2> "+tmp_dir+"/logminimap.txt" 
           # + " | awk 'BEGIN{FS=OFS=\"\t\"} {if($1 ~ /^@/) print $0; else {for(i=1; i<=NF; i++) {if(i==10) $i=\"*\";} print $0}}' > " + reads_on_asm + " 2> "+tmp_dir+"/logminimap.txt"
        print(" - Running minimap with command line:\n     " , command , "\n   The log of minimap2 can be found at "+tmp_dir+"/logminimap.txt")
        #write in the log file the time at which the alignment starts
        f = open(logFile, "a")
        f.write(" - Aligning the reads on the assembly\n")
        f.close()
        res_minimap = os.system(command)
        if res_minimap != 0 :
            print("ERROR: minimap2 failed. Was trying to run: " + command)
            print("ERROR: minimap2 could not run properly, check "+tmp_dir+"/logminimap.txt")
            sys.exit(1)

        #write in log file that alignment went smoothly
        f = open(logFile, "a")
        f.write("\nSTAGE 2: Alignment computed, minimap2 exited successfully\n")
        f.close()
    else:
        print(" - Already aligned reads found from previous run")

    print("\n===== STAGE 3: Calling variants   [", datetime.datetime.now() ,"]\n")
    sys.stdout.flush()

    #write in the log file the time at which the variant calling starts

    error_rate_file = tmp_dir + "/error_rate.txt"
    flag_debug = "0"
    if args.debug:
        flag_debug = "1"
    if continue_from_previous_run and os.path.exists(tmp_dir + "/variants.col") and os.path.exists(error_rate_file):
        print(" - Already called variants found from previous run")
    else:
        continue_from_previous_run = False
        command = path_call_variants + " " + new_assembly + " " + readsFile + " " + reads_on_asm + " " + str(nb_threads) + " " + tmp_dir + " " + error_rate_file + " " \
            + flag_debug + " " + tmp_dir + "/variants.col " + tmp_dir + "/variants.vcf"
        f = open(logFile, "a")
        f.write("\n==== STAGE 3: Calling variants   ["+str(datetime.datetime.now())+"]\n")
        f.write(command+"\n")
        f.close()
        # print(" - Calling variants with a basic pileup")
        print(" Running: ", command)
        res_call_variants = os.system(command)
        if res_call_variants != 0:
            print("ERROR: call_variants failed. Was trying to run: " + command)
            sys.exit(1)

        #write in the log file that variant calling went smoothly
        f = open(logFile, "a")
        f.write("STAGE 3: Variant calling computed, call_variants exited successfully. Variants are stored in "+tmp_dir+"/variants.vcf and "+tmp_dir+"/variants.col\n")
        f.close()

    #reading the error rate
    error_rate = 0.0
    with open(error_rate_file, 'r') as f:
        error_rate = float(f.readline())

    if error_rate > 0.15 :
        error_rate = 0.15 #more errors than this are probably heterozygous variants

    #write in the log file the error rate
    f = open(logFile, "a")
    f.write("STAGE 3: Error rate estimated from the alignment, error rate is "+str(error_rate))
    f.close()


    print("\n===== STAGE 4: Separating reads by haplotype of origin   [", datetime.datetime.now() ,"]\n")
    sys.stdout.flush()

    #estimate the ploidy of all the contigs if --haploid-coverage is used
    if haploid_coverage > 0 :
        if not continue_from_previous_run or not os.path.exists(tmp_dir + "/ploidy.txt"):
            continue_from_previous_run = False
            print(" - Estimating the ploidy of the contigs")
            command = path_determine_multiplicity + " " + new_assembly + " " + str(haploid_coverage) + " " + tmp_dir + "/ploidy.txt"
            print(" Running: ", command)
            res_estimate_ploidy = os.system(command)
            if res_estimate_ploidy != 0:
                print("ERROR: estimate_ploidy.py failed. Was trying to run: " + command)
                sys.exit(1)

            #write in the log file that ploidy estimation went smoothly
            f = open(logFile, "a")
            f.write("STAGE 4: Ploidy estimation computed, estimate_ploidy.py exited successfully. Ploidy is stored in "+tmp_dir+"/ploidy.txt")
            f.close()
    else:
        #create empty ploidy file
        f = open(tmp_dir + "/ploidy.txt", "w")
        f.close()

    #"Usage: ./separate_reads <columns> <num_threads> <error_rate> <DEBUG> <outfile> "
    command = path_separate_reads + " " + tmp_dir + "/variants.col " + str(nb_threads) + " " + str(error_rate) + " hs/tmp/ploidy.txt " + str(int(low_memory)) \
        + " " + str(rarest_strain_abundance) + " " + tmp_dir + "/reads_haplo.gro " + flag_debug
    #write in the log file the time at which the separation starts
    f = open(logFile, "a")
    f.write("\n==== STAGE 4: Separating reads by haplotype of origin   ["+str(datetime.datetime.now())+"]\n")
    f.write(command)
    f.close()

    if continue_from_previous_run and os.path.exists(tmp_dir + "/reads_haplo.gro") :
        print(" - Already separated reads found from previous run")
    else:
        continue_from_previous_run = False
        print(" - Separating reads by haplotype of origin")
        print(" Running: ", command)
        res_separate_reads = os.system(command)
        if res_separate_reads != 0:
            print("ERROR: separate_reads failed. Was trying to run: " + command)
            sys.exit(1)

        #write in the log file that read separation went smoothly
        f = open(logFile, "a")
        f.write("STAGE 4: Read separation computed, separate_reads exited successfully. Groups of reads are stored in "+tmp_dir+"/reads_haplo.gro. Explanation of the format\
                can be found in the doc/README.md, and a synthetic summary is in hairsplitter_summary.txt")

    print("\n===== STAGE 5: Creating all the new contigs   [", datetime.datetime.now() ,"]\n\n This can take time, as we need to polish every new contig using Racon")
    sys.stdout.flush()
    #"Usage: ./create_new_contigs <original_assembly> <reads_file> <error_rate> <split_file> <tmpfolder> <num_threads> <technology> <output_graph> <output_gaf> <MINIMAP> <RACON> <python> <debug>" 
    
    gaffile = tmp_dir + "/reads_on_new_contig.gaf"
    zipped_GFA = tmp_dir + "/zipped_assembly.gfa"
    polish_everything = "0"
    if args.polish_everything:
        polish_everything = "1"
    command = path_create_new_contigs + " " \
        + new_assembly + " " \
        + readsFile + " " \
        + str(error_rate) + " " \
        + tmp_dir + "/reads_haplo.gro " \
        + reads_on_asm + " " \
        + tmp_dir + " " \
        + str(nb_threads) + " " \
        + technology + " " \
        + zipped_GFA + " " \
        + gaffile +  " " \
        + polisher + " " \
        + polish_everything + " " \
        + path_to_minimap2 + " " \
        + args.path_to_racon + " " \
        + args.path_to_medaka + " " \
        + args.path_to_samtools + " " \
        + path_to_python + " " \
        + flag_debug
    print(" Running : ", command)
    #write in the log file the time at which the new contigs creation starts
    f = open(logFile, "a")
    f.write("\n==== STAGE 5: Creating all the new contigs   ["+str(datetime.datetime.now())+"]\n")
    f.write(command)
    f.close()

    if continue_from_previous_run and os.path.exists(zipped_GFA) :
        print(" - Already created new contigs found from previous run")
    else:
        continue_from_previous_run = False
        res_create_new_contigs = os.system(command)
        if res_create_new_contigs != 0:
            print("ERROR: create_new_contigs failed. Was trying to run: " + command)
            sys.exit(1)

    #write in the log file that new contigs were created
    f = open(logFile, "a")
    f.write("STAGE 6: New contigs created, create_new_contigs exited successfully. The new assembly graph is stored in "+zipped_GFA+" and the alignments of the reads\
            on the new contigs are stored in "+gaffile)
    f.close()

    print("\n===== STAGE 6: Untangling (~scaffolding) the new assembly graph to improve contiguity   [", datetime.datetime.now() ,"]\n")
    sys.stdout.flush()

    simply = ""
    if args.dont_simplify :
        simply = " --dont_merge"

    outfile = args.output.rstrip('/') + "/hairsplitter_final_assembly.gfa"


    command = path_graphunzip + " unzip -R -l " + gaffile + " -g " + zipped_GFA + simply + " -o " + outfile + " -r " + readsFile \
          + " 2>"+tmp_dir+"/logGraphUnzip.txt >"+tmp_dir+"/trash.txt"
    #write in the log file the time at which the untangling starts
    f = open(logFile, "a")
    f.write("\n==== STAGE 6: Untangling (~scaffolding) the new assembly graph to improve contiguity   ["+str(datetime.datetime.now())+"]\n")
    f.write(command)
    f.close()
    print( " - Running GraphUnzip with command line:\n     ", command, "\n   The log of GraphUnzip is written on ",tmp_dir+"/logGraphUnzip.txt\n")

    if continue_from_previous_run and os.path.exists(outfile) :
        print(" - Already untangled assembly found from previous run")
    else:
        continue_from_previous_run = False
        resultGU = os.system(command)
        if resultGU != 0 :
            print( "ERROR: GraphUnzip failed. Please check the output of GraphUnzip in "+tmp_dir+"/logGraphUnzip.txt" )
            sys.exit(1)
    
    #write in the log file that untangling went smoothly
    f = open(logFile, "a")
    f.write("STAGE 7: Untangling computed, GraphUnzip exited successfully. The new assembly is stored in "+outfile+". To see how the contigs were merged, check out hairsplitter_summary.txt.")
    f.close()

    print( "\n *To see in more details what supercontigs were created with GraphUnzip, check the hairsplitter_summary.txt*\n")
    output_file = "output.txt"
    o = open(output_file, "a")
    o.write("\n\n *****Linking the created contigs***** \n\nLeft, the name of the produced supercontig. Right, the list of new contigs with a suffix -0, -1...indicating the copy of the contig, linked with _ \n\n")
    o.close()
    command = "cat output.txt "+args.output +"/supercontigs.txt > output2.txt 2> "+args.output+"/tmp/trash.txt"
    os.system(command)
    command =  "mv output2.txt "+args.output+"/hairsplitter_summary.txt && rm supercontigs.txt output.txt 2> "+args.output+"/tmp/trash.txt";
    os.system(command)

    #write in the log file that the summary file was created
    f = open(logFile, "a")
    f.write("STAGE 6: Summary file created, hairsplitter_summary.txt is stored in "+args.output)
    f.close()

    
    fasta_name = outfile[0:-4] + ".fasta"
    command = path_gfa2fa + " " + outfile + " > " + fasta_name
    res_gfa2fasta = os.system(command)
    if res_gfa2fasta != 0:
        print("ERROR: gfa2fa failed. Was trying to run: " + command)
        sys.exit(1)

    if clean_tmp :
        command = "rm -r " + reads_on_asm + " " + tmp_dir + "/variants.col " + tmp_dir + "/variants.vcf " + tmp_dir + "/reads_haplo.gro " + tmp_dir + "/ploidy.txt " + tmp_dir + "/reads.fasta " + tmp_dir + "/reads_on_new_contig.gaf " + tmp_dir + "/reads_on_new_contig.gro " + tmp_dir + "/reads_on_new_contig.sam " + tmp_dir + "/reads_on_new_contig.bam " + tmp_dir + "/reads_on_new_contig.bam.bai " 
        res_clean = os.system(command)
        if res_clean != 0:
            print("ERROR: Could not remove temporary files. Was trying to run: " + command)
            #sys.exit(1)

    print("\n===== HairSplitter finished! =====   [", datetime.datetime.now() ,"]\n")

    #write in the log file that hairsplitter finished
    f = open(logFile, "a")
    f.write("\n==== HairSplitter finished!   ["+str(datetime.datetime.now())+"]\n")
    f.close()

if __name__ == "__main__":
    main()


