#!/usr/bin/env python

"""
HairSplitter takes as input a collapsed assembly and a set of long reads and outputs an uncollapsed assembly.
This is the master file that calls all the other scripts.
Author: Roland Faure

"""

__author__ = "Roland Faure"
__license__ = "GPL3"
__version__ = "1.5.0"
__date__ = "2023-10-18"
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
    parser.add_argument("-f", "--fastq", help="Sequencing reads fastq (required)", required=True)
    parser.add_argument("-x", "--technology", help="{ont, pacbio, hifi} [ont]", default="ont")
    parser.add_argument("-p", "--polisher", help="{racon,medaka,auto} medaka is more accurate but much slower. auto uses medaka on genomes <100kb, racon elsewhise [auto]", default="auto")
    # parser.add_argument("-m", "--multiploid", help="Use this option if all haplotypes can be assumed to have the same coverage", action="store_true")
    parser.add_argument("-t", "--threads", help="Number of threads [1]", default=1)
    parser.add_argument("-o", "--output", help="Output directory", required=True)
    parser.add_argument("-s", "--dont_simplify", help="Don't rename the contigs and don't merge them", action="store_true")
    parser.add_argument("-P", "--polish-everything", help="Polish every contig with racon, even those where there is only one haplotype ", action="store_true")
    parser.add_argument("-F", "--force", help="Force overwrite of output folder if it exists", action="store_true")
    parser.add_argument("--path_to_minimap2", help="Path to the executable minimap2 [minimap2]", default="minimap2")
    parser.add_argument("--path_to_racon", help="Path to the executable racon [racon]", default="racon")
    parser.add_argument("--path_to_medaka", help="Path to the executable medaka [medaka]", default="medaka")
    parser.add_argument("--path_to_samtools", help="Path to samtools [samtools]", default="samtools")
    parser.add_argument("--path_to_python", help="Path to python [python]", default="python")
    parser.add_argument("-d", "--debug", help="Debug mode", action="store_true")

    parser.add_argument("-v", "--version", help="Print version and exit", action="store_true")

    return parser.parse_args()

def check_dependencies(tmp_dir, minimap2, racon, medaka, polisher, samtools, path_to_src, path_to_python):

    com = " --version > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
    mini_run = os.system(minimap2 + com)
    if mini_run != 0:
        print("ERROR: minimap2 could not run. Check the path to the executable. (command line tried by HairSplitter: "+minimap2+com+")")
        sys.exit(1)

    if polisher != "medaka" :
        racon_run = os.system(racon + com)
        if racon_run != 0:
            print("ERROR: racon could not run. Check the path to the executable. (command line tried by HairSplitter: "+racon+com+")\n",
                  "If you want to use medaka, force the use of medaka through the --polisher option.")
            sys.exit(1)

    if polisher != "racon" :
        medaka_run = os.system(medaka + com)
        if medaka_run != 0:
            print("ERROR: medaka could not run. Check the path to the executable. (command line tried by HairSplitter: "+medaka+com+")\n",
                  "You can force the use of Racon instead, through the --polisher option.")
            sys.exit(1)

    samtools_run = os.system(samtools + com)
    if samtools_run != 0:
        print("ERROR: samtools could not run. Check the path to the executable. (command line tried by HairSplitter: "+samtools+com+")")
        sys.exit(1)

    command = path_to_python + " --version > "+tmp_dir+"/dependancies_log.txt 2> "+tmp_dir+"/dependancies_log.txt"
    python_run = os.system(command)
    if python_run != 0:
        print("ERROR: python could not run. Check the path to the executable. (command line tried by HairSplitter: "+command+")")
        sys.exit(1)

    command = path_to_python + " " + path_to_src+"GraphUnzip/graphunzip.py --help > "+tmp_dir+"/trash.txt 2> "+tmp_dir+"/dependancies_log.txt"
    graphunzip_run = os.system(command)
    if graphunzip_run != 0:
        print("ERROR: graphunzip could not run, probably because of a missing python package (check "+tmp_dir+"/dependancies_log.txt).\nCommand line tried by HairSplitter: "+command)
        sys.exit(1)

def main():

    if sys.argv[1] == "-v" or sys.argv[1] == "--version":
        print("HairSplitter v"+__version__+" ("+__github__+"). Last update: "+__date__)
        sys.exit(0)

    args = parse_args()
    nb_threads = args.threads
    #path to src folder can be imputed from the first argument of the command line
    path_to_src = sys.argv[0].split("hairsplitter.py")[0]+"src/"
    path_to_minimap2 = args.path_to_minimap2
    readsFile = args.fastq
    tmp_dir = args.output.rstrip('/') + "/tmp"
    path_to_python = args.path_to_python

    polisher = args.polisher.lower()
    if polisher != "racon" and polisher != "medaka" and polisher != "auto":
        print("ERROR: polisher must be either racon, medaka or auto")
        sys.exit(1)

    reads_on_asm = tmp_dir + "/reads_on_asm.sam"

    #print the command line used to run HairSplitter
    print(" ".join(sys.argv))
    print("HairSplitter v"+__version__+" ("+__github__+"). Last update: "+__date__)
    if args.version:
        sys.exit(0)

    #check if all the files and dependencies are here
    # check if output folder exists
    if os.path.exists(args.output) and not args.force:
        print("ERROR: output folder already exists. Use -F to overwrite.")
        sys.exit(1)
    elif not os.path.exists(args.output) :
        # create output folder
        os.mkdir(args.output)

    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    # check if input files exist
    if not os.path.exists(args.assembly):
        print("ERROR: not found assembly (" + args.assembly + ")")
        sys.exit(1)
    if not os.path.exists(args.fastq):
        print("ERROR: not found fastq file (" + args.fastq + ")")
        sys.exit(1)

    #check the dependencies
    check_dependencies(tmp_dir, args.path_to_minimap2, args.path_to_racon, args.path_to_medaka, args.polisher, args.path_to_samtools, path_to_src, path_to_python)

    # run the pipeline
    print("\n\t******************\n\t*                *\n\t*  Hairsplitter  *\n\t*    Welcome!    *\n\t*                *\n\t******************\n")

    logFile = args.output.rstrip('/') + "/hairsplitter.log"

    # 0. Convert the assembly to gfa if needed
    if args.assembly[-3:] == "gfa":
        gfaAssembly = args.assembly
    elif args.assembly[-5:] == "fasta" or args.assembly[-2:] == "fa":
        gfaAssembly = tmp_dir + "/assembly.gfa"
        command = path_to_src + "build/fa2gfa " + args.assembly + " > " + gfaAssembly
        res_fasta2gfa = os.system(command)
        if res_fasta2gfa != 0:
            print("ERROR: Conversion from fasta to gfa failed while running the command:\n" + command)
            sys.exit(1)
    else:
        print("ERROR: Assembly file must be in GFA or FASTA format. File extension not recognized.")
        sys.exit(1)

    # 1. Clean the assembly
    print("\n===== STAGE 1: Cleaning graph of small contigs that are unconnected parts of haplotypes   [", datetime.datetime.now() ,"]\n\n")
    print(" When the assemblers manage to locally phase the haplotypes, they sometimes assemble the alternative haplotype as a separate contig, unconnected in "\
                    "the gfa graph. This affects negatively the performance of Hairsplitter. Let's delete these contigs\n")
    
    #write in the log file where to look in case of error
    f = open(logFile, "w")
    f.write("==== STAGE 1: Cleaning graph of small contigs that are unconnected parts of haplotypes   ["+str(datetime.datetime.now())+"]\nIn case of errors\
            , check the log file of minimap2 at "+tmp_dir+"/tmp/minimap2.log")
    f.close()
    
    print(" - Mapping the assembly against itself")
    new_assembly = tmp_dir + "/cleaned_assembly.gfa"
    command = path_to_src + "build/clean_graph " + gfaAssembly + " " + new_assembly + " " + args.output.rstrip('/') +" " + logFile + " " \
        + str(nb_threads) + " " + args.path_to_minimap2
    print(" Running: ", command)
    res_clean = os.system(command)

    if res_clean != 0:
        print("ERROR: Cleaning the assembly failed. Was trying to run: " + command)
        sys.exit(1)

    print(" - Eliminated small unconnected contigs that align on other contigs")

    # 2. Map the reads on the assembly
    print("\n===== STAGE 2: Aligning reads on the reference   [", datetime.datetime.now() ,"]\n")

    #write in the log file the time at which the alignment starts
    f = open(logFile, "a")
    f.write("\n==== STAGE 2: Aligning reads on the reference   ["+str(datetime.datetime.now())+"]\n")
    f.close()

    # 2.1. Cut the contigs in chunks of 300000bp to avoid memory issues
    print(" - Cutting the contigs in chunks of 300000bp to avoid memory issues")
    command = "python " + path_to_src + "cut_gfa.py -a " + new_assembly + " -l 300000 -o " + tmp_dir + "/cut_assembly.gfa"
    res_cut_gfa = os.system(command)
    if res_cut_gfa != 0 :
        print("ERROR: cut_gfa.py failed. Was trying to run: " + command)
        sys.exit(1)
    new_assembly = tmp_dir + "/cut_assembly.gfa"

    # 2.2 Convert the assembly in fasta format
    print(" - Converting the assembly in fasta format")
    fastaAsm = tmp_dir + "/cleaned_assembly.fasta"
    command = path_to_src + "build/gfa2fa " + new_assembly + " > " + fastaAsm
    res_gfa2fasta = os.system(command)
    if res_gfa2fasta != 0 :
        print("ERROR: gfa2fa failed UUE. Was trying to run: " + command)
        sys.exit(1)
    
    # 2.3 Align the reads on the assembly
    print(" - Aligning the reads on the assembly")
    techno_flag = ""
    technology = args.technology.lower()
    if technology == "pacbio" or technology == "pb":
        techno_flag = "-x map-pb"
    elif technology == "hifi" :
        techno_flag = "-x map-hifi"
    else :
        techno_flag = "-x map-ont"

    command = path_to_minimap2 + " " + fastaAsm + " " + readsFile + " " + techno_flag + " -a -N 1 -t "+ str(nb_threads) +" > " + reads_on_asm + " 2> "+tmp_dir+"/logminimap.txt";
    print(" - Running minimap with command line:\n     " , command , "\n   The log of minimap2 can be found at "+tmp_dir+"/logminimap.txt")
    res_minimap = os.system(command)
    if res_minimap != 0 :
        print("ERROR: minimap2 failed. Was trying to run: " + command)
        print("ERROR: minimap2 could not run properly, check "+tmp_dir+"/logminimap.txt")
        sys.exit(1)

    #write in log file that alignment went smoothly
    f = open(logFile, "a")
    f.write("\nSTAGE 2: Alignment computed, minimap2 exited successfully")
    f.close()

    print("\n===== STAGE 3: Calling variants   [", datetime.datetime.now() ,"]\n")

    #write in the log file the time at which the variant calling starts
    f = open(logFile, "a")
    f.write("\n==== STAGE 3: Calling variants   ["+str(datetime.datetime.now())+"]\n")
    f.close()

    error_rate_file = tmp_dir + "/error_rate.txt"
    flag_debug = "0"
    if args.debug:
        flag_debug = "1"
    command = path_to_src + "build/call_variants " + new_assembly + " " + readsFile + " " + reads_on_asm + " " + str(nb_threads) + " " + tmp_dir + " " + error_rate_file + " " \
        + flag_debug + " " + tmp_dir + "/variants.col " + tmp_dir + "/variants.vcf"
    # print(" - Calling variants with a basic pileup")
    print(" Running: ", command)
    res_call_variants = os.system(command)
    if res_call_variants != 0:
        print("ERROR: call_variants failed. Was trying to run: " + command)
        sys.exit(1)

    #write in the log file that variant calling went smoothly
    f = open(logFile, "a")
    f.write("STAGE 3: Variant calling computed, call_variants exited successfully. Variants are stored in "+tmp_dir+"/variants.vcf and "+tmp_dir+"/variants.col")
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

    # print("\n===== STAGE 4: Filtering variants   [", datetime.datetime.now() ,"]\n")

    # #write in the log file the time at which the filtering starts
    # f = open(logFile, "a")
    # f.write("\n==== STAGE 4: Filtering variants   ["+str(datetime.datetime.now())+"]\n")
    # f.close()

    # command = path_to_src + "build/filter_variants " + tmp_dir + "/variants.col " + str(error_rate) + " " + str(nb_threads) + " " + flag_debug \
    #     + " " + tmp_dir + "/filtered_variants.col " + tmp_dir + "/variants.vcf " + tmp_dir + "/variants_filtered.vcf"
    # print(" - Filtering variants")
    # print(" Running: ", command)
    # res_filter_variants = os.system(command)
    # if res_filter_variants != 0:
    #     print("ERROR: filter_variants failed. Was trying to run: " + command)
    #     sys.exit(1)

    # #write in the log file that variant filtering went smoothly
    # f = open(logFile, "a")
    # f.write("STAGE 4: Variant filtering computed, filter_variants exited successfully. Filtered variants are stored in "+tmp_dir+"/variants_filtered.vcf \
    #         and "+tmp_dir+"/filtered_variants.col")
    # f.close()

    print("\n===== STAGE 4: Separating reads by haplotype of origin   [", datetime.datetime.now() ,"]\n")

    #write in the log file the time at which the separation starts
    f = open(logFile, "a")
    f.write("\n==== STAGE 4: Separating reads by haplotype of origin   ["+str(datetime.datetime.now())+"]\n")
    f.close()

    #"Usage: ./separate_reads <columns> <num_threads> <error_rate> <DEBUG> <outfile> "
    command = path_to_src + "build/separate_reads " + tmp_dir + "/variants.col " + str(nb_threads) + " " + str(error_rate) + " " + flag_debug \
        + " " + tmp_dir + "/reads_haplo.gro"
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
    #"Usage: ./create_new_contigs <original_assembly> <reads_file> <error_rate> <split_file> <tmpfolder> <num_threads> <technology> <output_graph> <output_gaf> <MINIMAP> <RACON> <python> <debug>" 
    
    #write in the log file the time at which the new contigs creation starts
    f = open(logFile, "a")
    f.write("\n==== STAGE 5: Creating all the new contigs   ["+str(datetime.datetime.now())+"]\n")
    f.close()

    if polisher == "auto" :
        if os.path.getsize(args.assembly) < 1000000 and args.technology == "ont" :
            polisher = "medaka"
        else :
            polisher = "racon"
    
    gaffile = tmp_dir + "/reads_on_new_contig.gaf"
    zipped_GFA = tmp_dir + "/zipped_assembly.gfa"
    polish_everything = "0"
    if args.polish_everything:
        polish_everything = "1"
    command = path_to_src + "build/create_new_contigs " \
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

    #write in the log file the time at which the untangling starts
    f = open(logFile, "a")
    f.write("\n==== STAGE 6: Untangling (~scaffolding) the new assembly graph to improve contiguity   ["+str(datetime.datetime.now())+"]\n")
    f.close()

    simply = ""
    if args.dont_simplify :
        simply = " --dont_merge -r"

    outfile = args.output.rstrip('/') + "/hairsplitter_final_assembly.gfa"

    meta = " --meta"
    # if args.multiploid :
    #     meta = ""

    command = "python " + path_to_src + "GraphUnzip/graphunzip.py unzip -l " + gaffile + " -g " + zipped_GFA + simply + " -o " + outfile + " 2>"+tmp_dir+"/logGraphUnzip.txt >"+tmp_dir+"/trash.txt";
    print( " - Running GraphUnzip with command line:\n     ", command, "\n   The log of GraphUnzip is written on ",tmp_dir+"/logGraphUnzip.txt\n")
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
    command = path_to_src + "build/gfa2fa " + outfile + " > " + fasta_name
    res_gfa2fasta = os.system(command)
    if res_gfa2fasta != 0:
        print("ERROR: gfa2fa failed. Was trying to run: " + command)
        sys.exit(1)

    print("\n===== HairSplitter finished! =====   [", datetime.datetime.now() ,"]\n")

    #write in the log file that hairsplitter finished
    f = open(logFile, "a")
    f.write("\n==== HairSplitter finished!   ["+str(datetime.datetime.now())+"]\n")
    f.close()

if __name__ == "__main__":
    main()


