import pickle #reading and writing files
import os #to run command lines from python (for blast especially)
from Bio.Blast.Applications import NcbiblastnCommandline #to use blast
from Bio.Blast import NCBIXML
#import subprocess #to run command lines from python


# names containing underscores do not work well. Here is a function to get the '_' out of the gfa
def take_underscores_out_of_gfa_names(fileIn, fileOut) :
    
    fo = open(fileOut, 'w')
    fi = open(fileIn, 'r')
    
    for line in fi :
        ls = line.split('\t')
        if 'S' in ls[0] :
            
            ls[1] = ''.join(ls[1].split('_'))
        fo.write('\t'.join(ls))


def cut_chromosomes(solutionFile, chunks = 2000):
    
    os.system('mkdir cut')
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
                q = open('cut/'+line.strip('>').strip('\n')+'_cut.fasta', 'w')
                print('Cutting ', line.strip('>').strip('\n'))
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
        return queryfiles

def assign_a_chromosome_to_each_contig(assemblyFile, queryfiles, fileOut = 'assign.pickle', chunks = 2000):
    
    #here the reference will be the assemblyFile
    os.system('mkdir tmp')
    os.system('makeblastdb -dbtype nucl -parse_seqids -in ' + assemblyFile + ' -out tmp/dtb')
    print('Done building the database')
    
    assign = {}
    aFile = open(assemblyFile, 'r')
    
    for line in aFile :
        if '>' in line[0] :
            assign[line.strip('>').strip('\n')] = []
    
    print(assign.keys())
    
    #queryfiles = ['2_cut.fasta', '12_cut.fasta']
    
    #map the pieces of chromosome on the assembly file
    
    for query in queryfiles :
        
        contigsHere = []
        print('Looking at chromosome : ', query)
        blastn_cline = NcbiblastnCommandline(query='cut/'+query, db='tmp/dtb', evalue=0.000001, outfmt=5, out='tmp/blast_sol'+query.strip('.fasta') + '.xml', task="megablast", num_alignments=5, qcov_hsp_perc = 90, num_threads = 3)
        out, err = blastn_cline()
        print('Blast done')
        
        # result_handle = open('tmp/blast_sol'+query.strip('.fasta') + '.xml')
        # blast_records = NCBIXML.parse(result_handle)
        # 
        # nowhere = 0
        # total = 0
        # 
        # for blast_record in blast_records :
        #     for alignment in blast_record.alignments:
        #         nbhits = 0
        #         for hsp in alignment.hsps:
        #             #print("identities ", hsp.identities)
        #             if hsp.identities > 0.9*chunks and nbhits == 0:
        #                 nbhits += 1
        #                 #if query.rstrip('_cut.fasta') not in assign[alignment.title[:-19]] :
        #                 assign[alignment.title[:-19]] += [query.rstrip('_cut.fasta')]
        #                 print(alignment.title[:-19], " matches to ", query.rstrip('_cut.fasta'))
        #         
        #         if nbhits == 0 :
        #             nowhere += 1
        #         total += 1

        
        #print('In assign_a_chromosome_to_each_contig, haven\'t managed to map', nowhere/total*100, ' % of the chunks')
    
    fo = open(fileOut, 'wb')
    pickle.dump(assign, fo)
    
    #os.system('rm -r tmp')
    
def manual_check(queryfiles, chunks, fileOut):
    
    assign = {}
    for query in queryfiles :
        
        result_handle = open('tmp/blast_sol'+query.strip('.fasta') + '.xml')
        blast_records = NCBIXML.parse(result_handle)
        
        nowhere = 0
        total = 0
        
        for blast_record in blast_records :
            for alignment in blast_record.alignments:
                nbhits = 0
                for hsp in alignment.hsps:
                    #print("identities ", hsp.identities)
                    if hsp.identities > 0.9*chunks and nbhits == 0:
                        nbhits += 1
                        #if query.rstrip('_cut.fasta') not in assign[alignment.title[:-19]] :
                        if alignment.title[:-19] not in assign :
                            assign[alignment.title[:-19]] = []
                        assign[alignment.title[:-19]] += [query.rstrip('_cut.fasta')]
                        # print(alignment.title[:-19], " matches to ", query.rstrip('_cut.fasta'))
                
                if nbhits == 0 :
                    nowhere += 1
                total += 1
    
    todelete = []
    for i in assign.keys() :
        if len(assign[i]) == 1 :
            todelete += [i]
            
    for i in todelete:
        del assign[i]
    
    fo = open(fileOut, 'wb')
    pickle.dump(assign, fo)
    #print(assign)


def check_phasing(assigned, fastaFile, fileOut) : # the contigs of the fasta file should be merged with '_' between the names
    
    fo = open(fileOut, 'w')
    
    ass = {}
    for i in assigned.keys():
        
        ass[i] = []
        for c in set(assigned[i]) :
            ass[i] += [(c, assigned[i].count(c))]
    
    fi = open(fastaFile, 'r')
    
    phasingErrors = 0
    for line in fi :
        
        if '>' in line[0] :
            listOfContigs = line.strip('>').strip('\n').split('-')
            listOfContigs = [''.join(i.split('_')[1:]) for i in listOfContigs]
            #print(listOfContigs)
            
            chromosomes = set()
            allchromosomes = set()
            first = True            
            for contig in listOfContigs :
                
                
                if contig not in assigned :
                    #print('There is a problem in the way the contigs are named, possibly if they contain _ or - in their names')
                    r = 0
                else :
                    if assigned[contig] != [] :
                        
                        if first :
                            first = False
                            chromosomes = set(assigned[contig])
                            allchromosomes = set(assigned[contig])
                            print('\nIn contig ', line)
                            fo.write('\nIn contig '+ line)
                        elif len(chromosomes) != 0 :
                            chromosomes = chromosomes.intersection(set(assigned[contig]))
                            allchromosomes = allchromosomes.union(set(assigned[contig]))
                            if len(chromosomes) == 0 :
                                #print ('\nPhasing error detected at contig ', listOfContigs)
                                phasingErrors += 1
                        print(contig, ' ', ass[contig])
                        fo.write(contig+ ' '+str( ass[contig])+'\n')
            
            #print('Contig ', listOfContigs, ' is in chromosome ', chromosomes)
    
    phasingSuccesses = 0
    for i in assigned.keys():
        if len(assigned[i]) > 0 :
            phasingSuccesses += 1
    
    print('\n\nTo summarize ', phasingSuccesses, ' contigs contain no phasing errors and ', phasingErrors, ' contain phasing errors')
        
#first cut each chromosome in chunks
# queryfiles = cut_chromosomes('potato/solution_article/RH89-039-16_potato_genome_assembly.v3.fa', chunks = 1000)
#print(queryfiles) #give the value outputted here to query files in the future
# print('Done cutting')

#then assign to each contig of the original assembly which chromosome(s) it belongs to
queryfiles = ['chr11_1_cut.fasta', 'chr11_2_cut.fasta', 'chr12_1_cut.fasta', 'chr12_2_cut.fasta']
#manual_check(['chr10_1_cut.fasta', 'chr10_2_cut.fasta'], 10000, 'tmp/chr10a-chr10b.assign.pickle')
assign_a_chromosome_to_each_contig('potato/stuberosum.hifiasm_l0.p_utg.fasta', queryfiles, 'potato/assign.pickle', chunks = 10000)

#if you have the queryfiles and the assign.pickle already you can skip the two above functions.

#check if the new supercontigs are composed of contigs coming from only one chromosome
# 
# f = open('tmp/chr10a-chr10b.assign.pickle', 'rb')
# assigned = pickle.load(f)
# 
# print(assigned)
#print(assigned['utg000030l'])
# assigned['881'],
# assigned['854'],
# assigned['776'])
# check_phasing(assigned, 'potato/stuberosum.hifiasm_l0.p_utg.graphunzip_hic_accept0.40_reject0.10.merged.fasta', 'potato/manual_check_chr10a-chr10b.txt')

