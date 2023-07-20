#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 17:19:55 2020

File dedicated to the analysis of the HiC coverage of contigs
"""

import matplotlib.pyplot as plt
import numpy as np
import basic_functions as bf


def determine_HiC_coverage(
    hiccontactsfile, info_contig, fragment_list, header = True
):  # returns the number of HiC contacts per basepair

    with open(hiccontactsfile) as f:

        coverage = {}
        for c in info_contig :
            coverage[c[0]] = 0

        for line in f:

            line = line.strip("\n")
            line = line.split("\t")

            if not header:  # because the first line is a header
                contact = [int(line[0]), int(line[1]), int(line[2])]

                contig1 = fragment_list[contact[0]][0]
                contig2 = fragment_list[contact[1]][0]
                coverage[contig1] += contact[2]
                coverage[contig2] += contact[2]
                
            header = False

    # xrange = [i for i in range(len(coverage))]
    # y = []
    # for i in coverage :
    #     y.append(coverage[i])
        
    #plt.xlim([0,20000])
    #plt.hist(y, bins = 4000)
    # plt.ylim([0,1])

    return coverage


def determine_unconnected_contigs(hiccontactsfile, fragmentList):

    # first we are going to look at which contigs have no HiC contacts
    with open(hiccontactsfile) as f:

        does_a_contig_interact_with_this_one = [
            False for i in range(fragmentList[-1][0] + 1)
        ]
        for line in f:

            line = line.strip("\n")
            line = line.split("\t")

            if line[0] != "487796":  # because the first line is a header
                contact = [int(line[0]), int(line[1]), int(line[2])]

                contig1 = fragmentList[contact[0]][0]
                contig2 = fragmentList[contact[1]][0]

                does_a_contig_interact_with_this_one[contig1] = True
                does_a_contig_interact_with_this_one[contig2] = True

        contigUnconnected = []
        for i in range(len(does_a_contig_interact_with_this_one)):
            if does_a_contig_interact_with_this_one[i] == False:
                contigUnconnected += [i]

        print(contigUnconnected)

def restrictionSitesInEachContigs(genomeFastaFile, restrictionSequence) :
    restrictionSites = {}
    
    with open(genomeFastaFile, 'r') as f :
        name = ''
        seq = ''
        fc = ''

        for line in f :
            if line[0] == '>' :
                restrictionSites[name] = seq.count(restrictionSequence)
                name = line.strip('>').strip('\n')
                seq = ''
            else :
                # if '>' in line :
                #     print('WHAT ?')
                seq += line.strip('\n')
        restrictionSites[name] = seq.count(restrictionSequence)
                
        restrictionSites.pop('')
        
    return restrictionSites


def check_if_there_are_restriction_fragments_in_this_contig(
    contig, restrictionSiteSequence, gfaFile
):

    with open(gfaFile) as f:

        for line in f:
            ls = line.split('\t')
            if 'S' in ls[0] and contig in ls[1] :
                return ls[2].count(restrictionSiteSequence)
            
    print("There is a problem with the input contig")
    return 0


def correlation_GCcontent_HiCcoverage(coverage, genomeFastaFile, unconnectedContigs):

    GCcontentOfContigs = [-1] * len(coverage)

    with open(genomeFastaFile) as f:

        step = 0
        for line in f:

            if not ">" in line:
                GCcontent = (line.count("G") + line.count("C")) / len(line)
                GCcontentOfContigs[int(step / 2)] = GCcontent

            step += 1

        extract = [
            GCcontentOfContigs[i]
            for i in range(len(GCcontentOfContigs))
            if i in unconnectedContigs
        ]
        extractCoverage = [
            coverage[i] for i in range(len(coverage)) if i in unconnectedContigs
        ]
        # plt.scatter(extract, extractCoverage)
        plt.scatter(GCcontentOfContigs, coverage)
        plt.xlabel("GC content")
        plt.ylabel("Coverage")
        plt.ylim([0, 1])

def correlation_coverage_restrictionSites(coverage, restrictionSites):
    
    keys = coverage.keys()
    x = []
    y = []
    for k in keys :
        x.append(restrictionSites[k])
        y.append(coverage[k])
        
    plt.scatter(x,y, alpha = 0.5)
    plt.xlabel('Number of restriction sites')
    plt.ylabel('Coverage')
    
fragmentsFile = "Arabidopsis/Arabidopsis_hybrid/HiCmapping/fragments_list.txt"
matrixFile = "Arabidopsis/Arabidopsis_hybrid/HiCmapping/abs_fragments_contacts_weighted.txt"
contigFile = "Arabidopsis/Arabidopsis_hybrid/HiCmapping/info_contigs.txt"
fastaFile = "Arabidopsis/Arabidopsis_hybrid/assembly.fasta"
gfaFile = "Arabidopsis/Arabidopsis_hybrid/assembly_graph.gfa"

print(check_if_there_are_restriction_fragments_in_this_contig('contig_6405', 'GATC', gfaFile))

# coverage = determine_HiC_coverage(matrixFile, bf.read_info_contig(contigFile), bf.read_fragment_list(fragmentsFile))
# restrictionSites = restrictionSitesInEachContigs(fastaFile, 'GATC')
# correlation_coverage_restrictionSites(coverage, restrictionSites)

# coverage = bf.import_from_csv('listsPython/HiCcoverage.csv')
# coverage = [x[0] for x in coverage]

# coverageWithoutTheZeros = [x for x in coverage if x != 0]
# coverageLog = [np.log10(x) for x in coverageWithoutTheZeros]

# plt.hist(coverageLog, bins = 160)
# plt.xlabel('Number of HiC contacts per bp (log scale)')
# plt.ylabel('Number of contigs')
# plt.xlim([0,3])

# correlation_GCcontent_HiCcoverage(coverage, 'data/Assembly.fasta', unconnectedcontigs)
print("Finished")

