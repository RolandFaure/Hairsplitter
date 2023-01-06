#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 13:05:13 2020

In this file, test functions to test our algorithm
"""

import matplotlib.pyplot as plt
import random
import numpy as np
from scipy.signal import argrelextrema #to detect local extrema of arrays
import time
import pickle #reading and writing files
import os #to run command lines from python (for blast especially)
from Bio.Blast.Applications import NcbiblastnCommandline #to use blast
from Bio.Blast import NCBIXML
import subprocess #to run command lines from python

import input_output as io
from input_output import load_gfa
from input_output import export_to_GFA
from input_output import load_interactionMatrix
from solve_ambiguities import solve_ambiguities
from solve_ambiguities import intensity_of_interactions
from solve_ambiguities import merge_adjacent_contigs

from evaluate_solution import score_output
from evaluate_solution import draw_distance_HiCcontacts_correlation
from evaluate_solution import heat_solution
from evaluate_solution import simulated_annealing
from loops import flatten_loop

from copy import deepcopy
import scipy.integrate as integrate
from scipy import sparse

def print_blast_coverage(cov) :
    
    for ch in cov :
        plt.plot(cov[ch])

def get_names(gfaFile) :
    with open(gfaFile) as f :
        names = []
        lengths = []
        index = 0
        for line in f :
            ls = line.split('\t')
            if 'S' in ls[0] :
                names.append(ls[1])
                lengths.append(len(ls[2]))
                index += 1
        return names, lengths
      
def blast_find_contig(contig, gfaFile, database) :
    seq = ''
    with open(gfaFile) as f :
        for line in f :
            ls = line.split('\t')
            if contig in ls[1] :
                seq = ls[2]
                fo = open('trash_me.fa', 'w')
                fo.write('>'+contig+'\n' +seq)
                fo.close()
                break
    
    cmd = 'blastn -query trash_me.fa -db '+ database +' -outfmt 5 -out ttrash_me.xml -task megablast'
    os.system(cmd)
    
    result_handle = open('ttrash_me.xml')
    blast_records = NCBIXML.parse(result_handle)
    
    for blast_record in blast_records :
        checked = False
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:

                if hsp.align_length == hsp.identities and hsp.identities == blast_record.query_length:
                    checked = True
                    print(contig, hsp.sbjct_start, ' - ', hsp.sbjct_start+hsp.align_length)
                if hsp.identities >  0.98*hsp.align_length and hsp.identities > 0.95 * blast_record.query_length and not checked:
                    print(contig, hsp.sbjct_start, ' - ', hsp.sbjct_start+hsp.align_length)
                    
    os.remove('ttrash_me.xml')
    os.remove('trash_me.fa')
    
#function that gives the order of the fasta contigs on the original sequence
def blast_solution(solutionFile, fastaAssemblyDtb, outdirectory, chunks = 1000):
    
    with open(solutionFile) as sol :
        
        queryfiles = []
        q = open('rien', 'w')
        seq = ''
        bit = 0

        for line in sol :
            if '>' in line :
                seq = ''
                bit = 0
                q.close()
                queryfiles += [line.strip('>').strip('\n')+'_cut.fasta']
                q = open(outdirectory +'/'+line.strip('>').strip('\n')+'_cut.fasta', 'w')
            else :
                
                char = 0
                l = list(line)
                
                while len(l)-char+len(seq) > chunks :
                    
                    while len(seq) < chunks :
                        seq += l[char]
                        char += 1
                        
                    #print('s ' + seq)
                    q.write('>'+str(bit)+'\n')
                    q.write(seq + '\n')
                    bit += 1
                    seq = ''
                
                seq += "".join(list(line)[char:])
                
        os.remove('rien')
    
    contigsInChromosomes = []
    fo = open(outdirectory+'/sols.chr', 'w')
    for query in queryfiles :
        
        contigsHere = []
        print('Looking at chromosome : ', query)
        blastn_cline = NcbiblastnCommandline(query=outdirectory +'/'+query, db=fastaAssemblyDtb, evalue=0.0000001, outfmt=5, out=outdirectory+'/blast_sol'+query.strip('.fasta') + '.xml', task="megablast")
        out, err = blastn_cline()
        
        result_handle = open(outdirectory+'/blast_sol'+query.strip('.fasta') + '.xml')
        blast_records = NCBIXML.parse(result_handle)
        
        index = 0
        last_contig = '-1'
        seq = ''
        for blast_record in blast_records :
            new_contig = []
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    #if hsp.expect < 0.00001 and hsp.identities > 5000 :
                    #if hsp.identities >  0.98*hsp.align_length :
                    #if hsp.identities >  0.98*hsp.align_length and hsp.identities > 0.9 * blast_record.query_length:
                    if hsp.align_length == hsp.identities and hsp.identities == blast_record.query_length:
                        #print(index, ' hit there: ', alignment.title)
                        new_contig += [alignment.title.split(' ')[-1]]
                       
            if len(new_contig) >= 1 :
                
                s = "/".join(new_contig)
                if s != last_contig :
                #if True :
                    last_contig = s
                    seq += s + ' - '
                    
                    if len(new_contig) == 1 :
                        contigsHere.append(s)
                
            else :
                if last_contig == '?' :
                    seq += '? - '
                last_contig = '?'
                
        contigsInChromosomes.append(contigsHere)
        
        fo.write('Looking at : '+ query + '\n')
        fo.write(seq + '\n\n\n\n')
    
    #looking at which contigs are specific to each chromosome
    spec1 = []
    for i in contigsInChromosomes[0] :
        if i not in contigsInChromosomes[1] :
            spec1 += [i]
            
    spec2 = []
    for i in contigsInChromosomes[1] :
        if i not in contigsInChromosomes[0] :
            spec2 += [i]
    
    file = open('Escherichia_Coli/1a1k/sol.pickle', 'wb')
    pickle.dump([spec1,spec2], file)
    
def check_with_solution(gfaFile) :
    
    file = open('Escherichia_Coli/1a1k/sol.pickle', 'rb')
    [spec1, spec2] = pickle.load(file)
    
    file.close()
    
    allcontigs = []
    gfa = open(gfaFile, 'r')
    for line in gfa :
        ls = line.split('\t')
        if 'S' in ls[0] :
            contigs = ls[1].split('_')
            allcontigs.append([i.split('-')[0] for i in contigs])
    
    c = 0
    for supercontig in allcontigs :
        s1 = False
        s2 = False
        for contig in supercontig :
            if contig in spec1 :
                s1 = True
            if contig in spec2 : 
                s2 = True
        
        if s1 or s2:
            c += 1
        
        if s1 and s2 :
            print('Look at here : ', supercontig)
    print('In total, there are ', c, ' haplotigs')
        
#blast_solution('Escherichia_Coli/1a1k/diploid1a1k.fasta', 'Escherichia_Coli/1a1k/dtb/dbk63', 'Escherichia_Coli/1a1k', 1000)
#check_with_solution('Escherichia_Coli/1a1k/unzipped_merged.gfa') 
    
#function to test if the unzipped gfa maps on the solution
def blast_check(fastaFile, database) :
    
    blastn_cline = NcbiblastnCommandline(query=fastaFile, db=database, evalue=0.0000000001, outfmt=5, out="check_blast.xml", task="megablast")
    out, err = blastn_cline()
    
    result_handle = open('check_blast.xml')
    blast_records = NCBIXML.parse(result_handle)
    
    index = 0
    
    fo = open('assign_contigs.txt', 'w')
   #true_align = {}
    
    for blast_record in blast_records :
        checked = False
        
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
            
                #if hsp.expect < 0.00000000000001 and hsp.identities > 5000 :
                #if hsp.identities >  0.98*hsp.align_length :
                if hsp.identities > 0.8 * blast_record.query_length:
                #if hsp.align_length == hsp.identities and hsp.identities == blast_record.query_length:
                    
                    #true_align[names[index]] = (hsp.sbjct_start, hsp.sbjct_start+hsp.align_length)
                    s = blast_record.query + '\t' + alignment.title +'\n'
                    print(blast_record.query, ' against ', alignment.title, ' and its stats are : ', hsp.identities,  hsp.align_length, blast_record.query_length, hsp.expect)
                    fo.write(s)
                    
                    checked = True
                    # coverage[-1].append((alignment.title, hsp.sbjct_start, hsp.sbjct_start+hsp.align_length))
                    # if len(names[index].split('_')) > 1 :
                    #     print('Goooood.... ', names[index])
                    # # print('db length : ', blast_record.database_length)
                    # print("****Alignment****")
                    # print("sequence:", alignment.title)
                    # print("end and beginning of subj : ", hsp.sbjct_start, hsp.sbjct_start+hsp.align_length)
                    # print('query length : ', blast_record.query_length)
                    # print("e value:", hsp.identities, hsp.align_length, hsp.positives)
                    # print(hsp.query[0:75] + "...")
                    # print(hsp.match[0:75] + "...")
                    # print(hsp.sbjct[0:75] + "...")
                # elif hsp.identities >  0.98*hsp.align_length and hsp.identities > 0.9 * blast_record.query_length: #still take results for the coverage
                #     nohitscov.append( (alignment.title, hsp.sbjct_start, hsp.sbjct_start+hsp.align_length) )           
                # if blast_record.query == '6531_6532' :
                #     print (hsp.align_length, hsp.identities)
        if not checked:
            if blast_record.query_length > 400 :
                print('OO...', blast_record.query, blast_record.query_length)
        
        index += 1
        
    #print(coverage)
    #now let's look at the coverage of parts of the genome
    
    # cov = {}
    # for j in database_len :
    #    cov[j] = [0 for j in range(int(database_len[j]/1000))] 
    
    # for hits in coverage :
    #     nbhits = len(hits)
        
    #     for seg in hits :
    #         for i in range(int(seg[1]/1000), int(seg[2]/1000)) :
    #             cov[seg[0]][i] += 1/nbhits
                
    # print_blast_coverage(cov)
    
#function that checks phasing using Hi-Fi
def check_phasing(assignFile, outFile, listOfSegments):
     
    f = open(assignFile, 'r')
    assign = {}
    for line in f :
        ls = line.strip('\n').split('\t')
        if ls[0] in assign :
            assign[ls[0]] += [ls[1]]
        else :
            assign[ls[0]] = [ls[1]]
            
    fo = open(outFile, 'w')
            
    for i in assign :
        i = list(set(i))
    
    for segment in listOfSegments :
        names = [i.split('-')[0] for i in segment.names[0].split('_')]
        
        bigcontig = [-1]
        
        for n in names :
            if n in assign :
                fo.write(n + str(assign[n])+'\n')
                print(n, ' ', assign[n])
                if bigcontig == [-1] :
                    bigcontig = assign[n]
                else :
                    newbig = []
                    for i in assign[n]:
                        if i in bigcontig :
                            newbig += [i]
                    bigcontig = newbig
            else :
                fo.write(n + ' :s\n')
                print(n + ' :s\n')
          
        if bigcontig == [] :
            fo.write(segment.full_name()+ ' is in contig : ' + str(bigcontig) + '\n\n')
            print('Contig ', segment.full_name(), ' is in contig : ', bigcontig, '\n\n')
        else :
            fo.write(segment.full_name()+ ' is in contig : ' + str(bigcontig) + '\n\n')
            print('Contig ' , segment.full_name(), '\n\n')
                    
    #print(assign)


bad_contigs = ['6', '70', '82', '87', '111', '119', '132', '144', '154', '166', '185', '195', '212', '213', '218', '229', '237']
long_contigs = ['39', '52', '53', '59', '60', '61', '62', '63', '64', '65', '66', '67', '70', '75', '77', '87', '92', '96', '103', '104', '115', '116', '118', '119', '120', '121', '125', '126', '132', '142', '144', '145', '152', '166', '184', '193', '194', '195', '201', '203', '211', '221', '224', '228', '230', '233', '234', '239']
true_align = {'0': (2215484, 2217866), '1': (2900783, 2902069), '10': (1271306, 1274458), '13': (1781175, 1785652), '14': (1562705, 1568145), '17': (2207382, 2209591), '18': (1657665, 1665888), '19': (5350691, 5351803), '20': (1558604, 1565501), '22': (2185577, 2187667), '24': (2913070, 2913783), '26': (2227469, 2228490), '27': (1181516, 1182446), '28': (1943491, 1944162), '29': (3482676, 3483764), '30': (4586843, 4587773), '32': (1181577, 1188070), '33': (1943552, 1944942), '34': (2695698, 2696295), '36': (2671982, 2672893), '39': (2708998, 2719781), '40': (2216113, 2222281), '41': (1558367, 1559962), '42': (1797968, 1799830), '45': (5020195, 5021273), '47': (2695658, 2696936), '52': (1909846, 2010374), '53': (1271045, 1332431), '59': (4977054, 5024971), '60': (5018080, 5056981), '61': (4736377, 4882915), '62': (228103, 456189), '63': (4162903, 4508476), '64': (3450409, 3481329), '65': (2192375, 2204914), '66': (2243929, 2455267), '67': (1551175, 1651610), '68': (5019591, 5020681), '69': (4737888, 4738978), '71': (808144, 815088), '74': (1658896, 1663978), '75': (3872396, 4033841), '76': (2166088, 2167266), '77': (900127, 917761), '79': (2241717, 2243866), '80': (2698988, 2701675), '83': (1451439, 1453532), '84': (2745893, 2747986), '85': (1186886, 1188114), '88': (2901821, 2903761), '89': (2675752, 2679085), '92': (2604236, 2753584), '93': (1182314, 1187389), '95': (2899203, 2900003), '96': (917031, 1164404), '97': (1643816, 1649645), '98': (1792624, 1794066), '99': (2167763, 2170546), '100': (2903782, 2905950), '101': (2236790, 2238630), '102': (3485709, 3487714), '103': (1775364, 1886586), '104': (1965879, 1978120), '105': (2908923, 2911839), '106': (2172913, 2175838), '107': (2212504, 2213053), '108': (3483156, 3483830), '109': (2680079, 2689518), '110': (3871228, 3872066), '114': (1800326, 1803010), '115': (4160687, 4287274), '116': (3448285, 3806831), '118': (4738503, 4831988), '120': (4834999, 4932471), '121': (5020206, 5246274), '122': (2236003, 2239097), '123': (1559140, 1562704), '124': (2185041, 2190744), '125': (2047988, 2161509), '126': (618996, 811477), '127': (1796231, 1796822), '128': (2165212, 2166087), '129': (2222375, 2225133), '130': (2185954, 2187689), '133': (3483449, 3484005), '134': (1944922, 1946137), '135': (2180026, 2181175), '136': (1632403, 1633552), '137': (1965670, 1968650), '139': (3484665, 3486165), '140': (2908952, 2912844), '141': (1552200, 1553759), '142': (1175577, 1186713), '143': (2242904, 2244463), '145': (4588494, 4735169), '146': (617058, 619707), '147': (2046050, 2048699), '148': (3449858, 3451265), '149': (4736928, 4738335), '150': (5018631, 5020038), '151': (5349402, 5350130), '152': (2745195, 2891454), '155': (5018118, 5018905), '156': (4832819, 4834095), '157': (4736415, 4737202), '158': (4162865, 4165033), '159': (3450371, 3451769), '160': (5019624, 5022048), '161': (3448867, 3451291), '162': (4161269, 4163693), '163': (3872128, 3872759), '164': (1209519, 1209842), '167': (1162672, 1165102), '168': (2710768, 2713198), '169': (2234074, 2238565), '170': (1780009, 1784500), '171': (4586479, 4587388), '172': (3870803, 3871712), '173': (2903310, 2904284), '174': (3483841, 3484568), '175': (2905390, 2906565), '176': (2679321, 2682186), '177': (896813, 900611), '178': (3489381, 3491533), '179': (2696060, 2696727), '180': (1777118, 1777780), '181': (3483955, 3484664), '182': (2238752, 2239085), '183': (2919693, 2926283), '184': (5498350, 5646538), '186': (3484455, 3485382), '187': (1968410, 1969538), '188': (229674, 232161), '189': (3484496, 3485084), '190': (1943938, 1944838), '191': (4586975, 4588349), '192': (3871299, 3872673), '193': (2736004, 2762972), '194': (5217161, 5282001), '196': (2161181, 2162275), '197': (2670999, 2671572), '198': (1571971, 1576285), '199': (2671441, 2672482), '201': (2921080, 3089480), '202': (2208436, 2212232), '203': (1382692, 1490989), '204': (1807964, 1816022), '205': (4587895, 4589535), '206': (1954102, 1955313), '207': (1941488, 1942645), '208': (2227930, 2229047), '209': (2173389, 2174419), '210': (2909399, 2910429), '211': (273419, 346150), '214': (2683272, 2690162), '215': (1196072, 1202738), '216': (2898465, 2905217), '217': (2670219, 2670681), '219': (2901379, 2903285), '220': (2229582, 2230733), '221': (1450213, 1517612), '222': (2744667, 2751673), '223': (2242304, 2243664), '224': (5246451, 5349903), '225': (1633092, 1637812), '226': (2179337, 2184057), '228': (1197159, 1207973), '230': (1575310, 1633091), '231': (2166520, 2167818), '232': (2902539, 2903837), '233': (617989, 888662), '234': (2046981, 2124246), '235': (3483239, 3483864), '236': (1209604, 1211549), '238': (2178873, 2184682), '239': (1633556, 1644546), '240': (1791655, 1797219), '241': (2211779, 2212424)}

gfaFile = 'Escherichia_Coli/assemblyGraph_bwise_x_k63.gfa'
fastaFile = 'data_A_Vaga_PacBio/Assembly.fasta'
# gfaFile = 'Escherichia_Coli/unzipped1*.gfa'
# fastaFile = 'Escherichia_Coli/unzipped1*.fasta'
#names, lengths = get_names(gfaFile)
#segments, name = load_gfa('data_A_Vaga_PacBio/unzipped_merged.gfa')

#blast_check(fastaFile, 'data_A_vaga_HiFi/dtb/dtb_HiFi')
#blast_find_contig('221', gfaFile, 'Escherichia_Coli/databases/diploid1a1b.dtb')
#check_phasing('assign_contigs.txt', 'test_assigns.txt', segments)

def scan_interactionMatrix(gfa, matrix, contigs):
    segments, names = io.load_gfa(gfa)
    interactions = io.load_interactionMatrix(matrix, segments, names)
    
    for i in contigs :
        print('scan_interactions : ', interactions[names[i[0]],names[i[1]]])

# gfaFile = 'data_A_vaga_HiFi/Flye/assemblyFlyeHiFi.gfa'
# matrix = 'data_A_vaga_HiFi/Flye/interactionMatrix.pickle'
# 
# contigs = [['edge_267', 'edge_78'], ['edge_268', 'edge_78']]
# scan_interactionMatrix(gfaFile, matrix, contigs)

def testRatios():
    file = open('ratio.txt', 'r')
    ra = []
    for l in file :
        ra += [float(l)]
        
    plt.hist(ra, bins = 40)
    plt.xlabel('Ratio between the two choices when the programm has to make a choice')
    plt.ylabel('Number of occurences')

def print_chromosomes(chromosomes):
    for c in chromosomes :
        s=''
        for contig in c :
            s += '-'+contig
        print(s[1:])
        
def buildFakeChromosomes(chromosomesLength = 10):
    letters = ['A','A','B','B','C','C','D','D']
    chromosomes = [[letters[j]+str(i) for i in range(chromosomesLength)] for j in range(4)]
    
    duplicates = 0 #number of conitgs that are going to be repeated within each chromosomes :
    for chromosome in chromosomes :
        for i in range(duplicates):
            contigCopied = random.randint(0, len(chromosome)-1)
            insertionPosition = random.randint(0, len(chromosome)-1)
            chromosome.insert(insertionPosition, chromosome[contigCopied])
            
    mutations = 5 #number of total mutations in the genome 
    for i in range(mutations):
        c = random.randint(0, len(chromosomes)-1)
        p = random.randint(0, len(chromosomes[c])-1)
        chromosomes[c][p] += '*'
            
    crossCopies = 2 #number of contigs that are going to be randomly written somewhere in the genome
    for i in range(crossCopies):
        c = random.randint(0, len(chromosomes)-1)
        p = random.randint(0, len(chromosomes[c])-1)
        contigCopies = chromosomes[c][p]
        
        c = random.randint(0, len(chromosomes)-1)
        p = random.randint(0, len(chromosomes[c])-1)
        chromosomes[c].insert(c, contigCopies)
        
    repeatedElement = 0 #one transposon is going to go everywhere
    for i in range(repeatedElement):  
        c = random.randint(0, len(chromosomes)-1)
        p = random.randint(0, len(chromosomes[c])-1)
        chromosomes[c].insert(p, 'T')
    
    return chromosomes

def exportFakeToGFA(chromosomes, file, lengthOfContig) :
    
    f = open(file, 'w')
    
    alreadyExported = []
    for c in chromosomes :
        for contig in c :
            if contig not in alreadyExported :
                f.write('S\t'+contig+'\tCT'+ ''.join(['A' for i in range(lengthOfContig-2)])+'CT' +'\n')
                alreadyExported += [contig]
      
    linksAlreadyExported = []
    for c in chromosomes :
        for contig in range(len(c)-1) :
            if [c[contig],c[contig+1]] not in linksAlreadyExported :
                linksAlreadyExported += [[c[contig],c[contig+1]]]
                f.write('L\t'+c[contig]+'\t+\t'+c[contig+1]+'\t+\t2M\n')
                
    return alreadyExported
      
def dist_law(distance) :
    if distance < 1000 :
        distance = 1000
    return 1000000/(distance**2)

def constructFakeInteractionMatrix(chromosomes, names, segments, lengthOfContigs = 10000):
    
    interactionMatrix = sparse.dok_matrix((len(names), len(names)))
    for c in chromosomes :
        for c1 in range(len(c)):
            for c2 in range(c1, len(c)) :
                con1 = c[c1]
                con2 = c[c2]
                if con1 != 'T' and con2 != 'T' : #because 'T' elements are transposons
                    intensity = integrate.quad(dist_law, (c2-c1-1)*lengthOfContigs, (c2-c1)*lengthOfContigs)[0]
                    interactionMatrix[names[con1],names[con2]] += intensity
                    interactionMatrix[names[con2],names[con1]] += intensity
                    segments[names[con1]].HiCcoverage += intensity
                    segments[names[con2]].HiCcoverage += intensity
    
    return interactionMatrix

#function that checks if ls1 is a sublist of ls2 : useful for checking the quality of the output
def sublist(ls1, ls2):
    
    sub = False
    
    s1 = ','.join(ls1)
    s2 = ','.join(ls2)
    sub = sub or (s1 in s2)
    
    s1 = ','.join(ls1[::-1]) #if it goes the other way, it is still imbricated
    s2 = ','.join(ls2)
    sub = sub or (s1 in s2)

    return sub

#function taking as arguments the solution of the problem and the output of the algorithm to see if the output is mistaken
def check_result(chromosomes, listOfSuperContigs, names) :

    #First check if all supercontigs of the output actually exist in the solution (i.e. contigs were not accidentally duplicated)
    for supercontig in listOfSuperContigs :
        
        found = False
        for c in chromosomes :
            if sublist(supercontig.names,c) :
                found = True
        if not found :
            print('Output : ')
            for i in listOfSuperContigs :
                i.print_complete()
            print('Actual chromosomes : ')
            for c in chromosomes :
                print(c)
            print('Contig found in output but not in chromosomes : ', supercontig.names)
            return False
        
    #Then check if all true links still exist (i.e. links were not accidentally deleted)

    expectedContacts = sparse.dok_matrix((len(names), len(names))) #this is the matrix of all the links that are in the solution and therefore are expected in the output
    for c in chromosomes :
        for contig in range(len(c)-1) :
            expectedContacts[names[c[contig]],names[c[contig+1]]] = 1
            expectedContacts[names[c[contig+1]],names[c[contig]]] = 1
 
        #First, look if the expected links are found between supercontigs
    for segment in listOfSuperContigs :
        for endOfSegment in range(2):
            for l, neighbor in enumerate(segment.links[endOfSegment]) :
                expectedContacts[ names[segment.names[-endOfSegment]] , names[neighbor.names[-segment.otherEndOfLinks[endOfSegment][l]]]] = 0

        #Then within supercontigs
    for segment in listOfSuperContigs :
        for j in range(len(segment.listOfContigs)-1) :
            expectedContacts[names[segment.names[j]] , names[segment.names[j+1]]] = 0
            expectedContacts[names[segment.names[j+1]] , names[segment.names[j]]] = 0
            
    for i in range(len(names)):
        for j in range(len(names)):
            if expectedContacts[i,j] == 1 :
                print('In the actual chromosomes, there exist a link between ', i, ' and ', j, ' that is not found in the output')
                print(names)
                for k in listOfSuperContigs :
                    k.print_complete()
                return False
    
    return True


def stats_on_solve_ambiguities(n = 100, lengthOfChromosomes = 10, steps = 10) :
    
    record = []
    for i in range(n):
        chromosomes = buildFakeChromosomes(lengthOfChromosomes)
        lengthOfContig = 10000
        exportFakeToGFA(chromosomes, 'tests/stats/test' + str(i)+'.gfa', lengthOfContig)
        bf.export_to_csv(chromosomes, 'tests/stats/test' + str(i)+'.chro')
        
        listOfSegments, names = load_gfa('tests/stats/test' + str(i)+'.gfa')

        interactionMatrix = constructFakeInteractionMatrix(chromosomes, names, lengthOfContig)

        listOfSegments = solve_ambiguities(listOfSegments, interactionMatrix, 0.2, 0.45 ,5) #rejectedThreshold<AcceptedThreshold

        record.append(check_result(chromosomes, listOfSegments, names))
        
        bf.export_to_GFA(listOfSegments, exportFile = 'tests/stats/test' + str(i)+'F.gfa')
        #draw_distance_HiCcontacts_correlation(listOfSuperContigs, links, [10000 for i in names], interactionMatrix)
        
    fileRecord = open('tests/stats/record.txt', 'w')
    for i in record :
        fileRecord.write(str(i)+'\n')
    print(int(record.count(False)*100/n), '% of incorrectly changed GFA')
    
def stats_on_thresholds(segments, names, interactionMatrix) :
    
    copiesNumber = {}
    for segment in segments :
        copiesNumber['_'.join(segment.names)] = 1

    ratios = []
    for segment in segments :
        
        for endOfSegment in range(2) :
            
            if len(segment.links[endOfSegment]) >= 2 : #then it means that there is a choice to be made at one end of the segment. Let's see how HiC contacts confirm those links
                    
                # comparison pairwise of the links, those that should be deleted are deleted
                    for n1 in range(len(segment.links[endOfSegment]) - 1):
                        n2 = n1 + 1
                        while n2 < len(segment.links[endOfSegment]):
                            
                            d = 2
                               
                            absoluteLinksStrength, linksStrength, neighborsOfNeighborsUsed = intensity_of_interactions(segment, [segment.links[endOfSegment][n1], segment.links[endOfSegment][n2]],\
                                                                                             [segment.otherEndOfLinks[endOfSegment][n1], segment.otherEndOfLinks[endOfSegment][n2]],\
                                                                                             segments, interactionMatrix, names, copiesNumber, depthOfCommonContigs = d)
                                                                                             
                            n2 += 1
                            
                            if len(linksStrength) == 2 and  np.max(linksStrength) > 2:
                                ratios += [ np.min(linksStrength)/np.max(linksStrength) ]
                                
    plt.hist(ratios)
    plt.xlabel('i(X)/i(Y) ratio')
    plt.ylabel('Number of ambiguities having this value')
    plt.show()
    
    return ratios
    
        
segments, names = load_gfa('data_A_vaga_HiFi/Flye/assemblyFlyeHiFi+.gfa')
interactionMatrix = load_interactionMatrix('data_A_vaga_HiFi/Flye/interactionMatrix.pickle', segments, names)
stats_on_thresholds(segments, names, interactionMatrix)

# t = time.time()
# # # chromosomes = ['A0-A1-A2-A3-A4-A5-A6-A7-A8-A9'.split('-'), 'A0-A1-A2-A3*-A4-A5-A6-A7-A8-A9'.split('-'),\
# # #                 'B0*-B1-B1-B2-B3-B4*-B5-B6-B7-B8-B9'.split('-'), 'B0*-B1-B2*-B3-B4-B5-B6-B7-B8-B9'.split('-')]

# # chromosomes = bf.import_from_csv('tests/fake.chro')
# chromosomes = buildFakeChromosomes(10)
# f = open('tests/fake.chro', 'wb')
# pickle.dump(chromosomes, f)

# lengthOfContig = 20
# exportFakeToGFA(chromosomes, 'tests/fake.gfa', lengthOfContig)
# listOfSegments, names = load_gfa('tests/fake.gfa')


# # interactionMatrix = constructFakeInteractionMatrix(chromosomes, names, listOfSegments, lengthOfContig)
# interactionMatrix = sparse.dok_matrix((len(names), len(names)))
# listOfSegments = solve_ambiguities(listOfSegments, interactionMatrix, names , 0.1, 0.2 ,3) #rejectedThreshold<AcceptedThreshold
# export_to_GFA(listOfSegments, 'tests/fake.gfa', 'tests/export.gfa', merge_adjacent_contigs = True)

# # export_to_GFA(listOfSegments, 'tests/fake.gfa', exportFile = 'tests/fakeF.gfa', merge_adjacent_contigs = False)

# # links, listOfSuperContigs, cn = simulated_annealing(originalLinks, names, interactionMatrix, [lengthOfContig for i in names], lambda x:1, 0.2, 0.45 ,5)
# # export_to_GFA(links, listOfSuperContigs, cn, originalLinks, names = names, exportFile = 'tests/fakeA.gfa')

# # print('Before the beginning of the process, the gfa energy is : ', \
# #       score_output([[i] for i in range(len(names))], originalLinks, [lengthOfContig for i in names], interactionMatrix, infinite_distance = 500000))

# #print('And the output is : ', check_result(chromosomes, listOfSegments, names))#, ', of energy ', score_output(listOfSuperContigs, links, [lengthOfContig for i in names], interactionMatrix, infinite_distance = 500000))

# #stats_on_solve_ambiguities(n=100)


# print('Finished in ', time.time()-t, ' seconds')






    