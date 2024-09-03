#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import re

def reverse_complement(seq) :
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    return "".join(complement[base] for base in seq[::-1])

#input: the graph (as the list of segments), the alignment of the reads (gaf_file), and the number of copies of each contig in the final assembly
#output: each subcontig assigned to a precise list of reads that make up the contig
def assign_reads_to_contigs(segments, gaf_file, copies):

    #go through the gaf file and associate reads to segments
    subcontig_to_segments = {}
    for s, segment in enumerate(segments) :
        for su, subcontig in enumerate(segment.get_namesOfContigs()) :
            if subcontig not in subcontig_to_segments :
                subcontig_to_segments[subcontig] = set()
            subcontig_to_segments[subcontig].add((s, su))
    
    #go through the gaf file and associate reads to segments
    with open(gaf_file, 'r') as gaf :
        for line in gaf :
            if line[0] == 'H' :
                continue
            line = line.strip().split('\t')
            read = line[0]
            path = line[5]
            contigs = re.split('[><]' , path)
            orientations = "".join(re.findall("[<>]", path))
            del contigs[0] #because the first element is always ''

            # print("Looking at path: ", path)

            # if read != "25B2_SRR19419643.41073":
            #     continue

            if len(contigs) > 0 :
                
                start_contig = 0
                for start_contig in range(len(contigs)) :
                    if contigs[start_contig] not in subcontig_to_segments:
                        continue

                    for c, candidate_segments in enumerate(subcontig_to_segments[contigs[start_contig]]) :
                        
                        unique_segment = [] #list of the subcontigs that are between two unique subcontigs (copies=1)
                        last_unique = -1
                        #check if the path matches the segment
                        good_path = True
                        index_in_contigs = start_contig

                        seg = segments[candidate_segments[0]]
                        index_in_names = candidate_segments[1]

                        names_seg = seg.get_namesOfContigs()
                        orientations_seg = seg.get_orientations()

                        inverse = 1
                        if orientations[start_contig] != "<>"[orientations_seg[index_in_names]] :
                            inverse = -1

                        #avoid checking the same path twice
                        if index_in_contigs == 0 or (index_in_names == 0 and inverse == 1) or (index_in_names == len(seg.get_namesOfContigs()) - 1 and inverse == -1) :
                        
                            if copies[contigs[index_in_contigs]] == 1 :
                                last_unique = index_in_contigs
                                unique_segment += [contigs[index_in_contigs]]
                            
                            while index_in_contigs < len(contigs) and index_in_names >= 0 and index_in_names < len(seg.get_namesOfContigs()) :
                                #see if we can continue the path
                                expected_orientation_in_seg = "<>".index(orientations[index_in_contigs])
                                if inverse == -1 :
                                    expected_orientation_in_seg = 1 - expected_orientation_in_seg

                                if contigs[index_in_contigs] != names_seg[index_in_names] or orientations_seg[index_in_names] !=  expected_orientation_in_seg :
                                    good_path = False
                                    break

                                else :
                                    if copies[contigs[index_in_contigs]] == 1 :
                                        if last_unique != -1 :
                                            unique_segment += contigs[last_unique+1:index_in_contigs+1]
                                        last_unique = index_in_contigs

                                    #move to the next contig
                                    index_in_contigs += 1
                                    index_in_names += inverse
                            
                            # print("was it a good path? ", good_path, " ", unique_segment)
                            if good_path :
                                # print("In segment: ", segments[candidate_segments[0]].get_namesOfContigs(), "\npath", path, "\nstarring ", unique_segment, "\nread ", read)
                                for contig in unique_segment :
                                    # print("adding read ", (contig, read), " to ", seg.names)
                                    seg.add_read(contig, read)

#input: the graph (as the list of segments), the alignment of the reads (gaf_file), and the number of copies of each contig in the final assembly and the fasta/q file and the gfa file
#output: repolished sequences stored in the subcontigs
def repolish_contigs(segments, gfa_file, gaf_file, fastq_file, copies, threads=1):

    #first assign all the reads to the subcontigs
    assign_reads_to_contigs(segments, gaf_file, copies)

    #index all of the reads with their position in the fasta/q file so that we can retrieve them later
    reads_position = {}
    line_number = 0
    last_record = -3
    fasta = False
    if fastq_file[-6:] == ".fasta" or fastq_file[-3:] == ".fa" :
        fasta = True
    with open(fastq_file, 'r') as fastq :
        line = fastq.readline()
        while line :
            if line[0] == '@' and line_number%4 == 0 and not fasta : #to avoid indexing quality scores
                read = line[1:].strip().split()[0]
                reads_position[read] = fastq.tell()
                fastq.readline()
                fastq.readline()
                line = fastq.readline()
                last_record = line_number
                line_number += 3
            elif line[0] == '>' and fasta:
                read = line.split()[0][1:].strip()
                reads_position[read] = fastq.tell()
                line = fastq.readline()
                line_number += 1
            else :
                line = fastq.readline()
                line_number += 1

    #index the contigs with their position in the gfa file so that we can retrieve them later
    contigs_position = {}
    previous_position = 0
    with open(gfa_file, 'r') as gfa :
        line = gfa.readline()
        while line :
            if line[0] == 'S' :
                line = line.strip().split('\t')
                contig = line[1]
                contigs_position[contig] = previous_position

            previous_position = gfa.tell()
            line = gfa.readline()

    #go through the segments and their subcontigs and repolish them using racon
    for segment in segments :
        seqs = segment.get_sequences()
        names = segment.get_namesOfContigs()
        reads = segment.get_reads()
        orientations = segment.get_orientations()

        for s, subcontig in enumerate(names) :

            # if subcontig != "edge_21_149088_182933_0_33845_0_33845@0_10000_0":
            # if "edge_18_310438_361584_0_51146_0_51146@0_10000_1" not in subcontig :
            #     print("continuuedj ", subcontig)
            #     continue

            # print("reads here: ",[len(i) for i in reads])
            # print(names[13], " ", reads[13])

            if segment.get_lengths()[s] < 100 :
                continue

            # print("Looking at subcontig ", subcontig, " ", s , " ", copies[subcontig], " ", len(reads[s]))
            if len(reads[s]) > 0 and copies[subcontig] > 1 : #if the contig is unique it should be already polished

                seq = None
                # print("Repolishing ", subcontig, " with ", len(reads[s]), " reads")

                #now repolish
                #begin by extracting the reads from the fastq file and write them to a temporary file
                f = open("tmp_reads.fa", 'w')
                with open(fastq_file, 'r') as fastq :
                    for read in reads[s] :
                        fastq.seek(reads_position[read])
                        f.write(">" + read + "\n")
                        f.write(fastq.readline())

                f.close()

                #find out the chunk of the contig left of the subcontig
                left = ""
                name_of_contig_left = names[s-1]
                with open(gfa_file, 'r') as gfa :
                    gfa.seek(contigs_position[name_of_contig_left])
                    ls = gfa.readline().strip().split('\t')
                    left = ls[2]
                    if orientations[s-1] == 0 :
                        left = reverse_complement(left)
                #write down left in a temporary file
                # print("left contig: ", name_of_contig_left)
                f = open("tmp_left.fa", 'w')
                f.write(">" + name_of_contig_left + "\n" + left + "\n")
                f.close()

                #find out the chunk of the contig right of the subcontig
                right = ""
                name_of_contig_right = names[s+1]
                with open(gfa_file, 'r') as gfa :
                    gfa.seek(contigs_position[name_of_contig_right])
                    ls = gfa.readline().strip().split('\t')
                    right = ls[2]
                    if orientations[s+1] == 0 :
                        right = reverse_complement(right)

                #write down right in a temporary file
                f = open("tmp_right.fa", 'w')
                f.write(">" + name_of_contig_right + "\n" + right + "\n")
                f.close()

                #first check that the reads align well on the contig - if not (e.g. structural variant), reassemble everythin
                contig_seq = ""
                contig_extended = ""
                with open(gfa_file, 'r') as gfa :
                    gfa.seek(contigs_position[subcontig])
                    ls = gfa.readline().strip().split('\t')
                    contig_seq = ls[2]
                    contig_extended = contig_seq
                    if orientations[s] == 0 : #if reverse complement
                        contig_extended = reverse_complement(contig_seq)
                    gfa.seek(0)
                    #if neighboring contigs are there let's take them too
                    if s > 0 and s < len(names)-1 :
                        gfa.seek(contigs_position[names[s-1]])
                        ls = gfa.readline().strip().split('\t')
                        neigh_seq = ls[2]
                        if orientations[s-1] == 0 : #if reverse complement
                            neigh_seq = reverse_complement(neigh_seq)
                        contig_extended = neigh_seq[-1000:] + contig_extended
                        gfa.seek(0)
                        gfa.seek(contigs_position[names[s+1]])
                        ls = gfa.readline().strip().split('\t')
                        neigh_seq = ls[2]
                        if orientations[s+1] == 0 : #if reverse complement
                            neigh_seq = reverse_complement(neigh_seq)
                        contig_extended = contig_extended + neigh_seq[:1000]
                f = open("tmp_complete_contig.fa", 'w')
                f.write(">" + subcontig + "_and_left_and_right" + "\n" + contig_extended + "\n")
                f.close()

                # align reads on the contig using minimap2
                command = "minimap2 -x map-pb -t " + str(threads) + " tmp_complete_contig.fa tmp_reads.fa > tmp_complete.paf 2> trash.txt"
                minimap = os.system(command)
                if minimap != 0 :
                    print("Error while running minimap2: " + command + "\n")
                    sys.exit(1)

                #check if the alignments (or at least one) are good
                no_struct_variants = False
                orientations_of_reads = {}
                with open("tmp_complete.paf", 'r') as paf :
                    for line in paf :
                        ls = line.strip().split('\t')
                        orientations_of_reads[ls[0]] = ls[4]
                        #make sure the read aligns on more or less the whole read
                        if int(ls[8])-int(ls[7]) > 0.9*int(ls[6]) \
                            and int(ls[7]) < 500 and int(ls[8]) > len(contig_seq)-500 \
                            and int(ls[3])-int(ls[2]) > 0.9*(int(ls[8])-int(ls[7])) and int(ls[3])-int(ls[2]) < 1.1*(int(ls[8])-int(ls[7])) :
                            no_struct_variants = True

                # print("no struct variants: ", no_struct_variants, " (", names[max(s-1, 0)], " ", names[min(s+1, len(names)-1)], ") ", s , " ", len(names)-1)

                if no_struct_variants or s == 0 or s == len(names)-1 :

                    print("polishing ", subcontig, " with ", len(reads[s]), " reads")
                    #output the contig to a temporary file
                    f = open("tmp_contig.fa", 'w')
                    f.write(">" + subcontig + "\n" + contig_seq + "\n")
                    f.close()

                    #now polish the contig with the reads using racon
                    command = "minimap2 -x map-pb -t " + str(threads) + " tmp_contig.fa tmp_reads.fa > tmp.paf 2> trash.txt"
                    minimap = os.system(command)
                    if minimap != 0 :
                        print("Error while running minimap2: " + command + "\n")
                        sys.exit(1)

                    command = "racon -t " + str(threads) + " tmp_reads.fa tmp.paf tmp_contig.fa > tmp_repolished.fa 2>trash.txt"
                    racon = os.system(command)
                    if racon != 0 :
                        # print("Error while running racon: " + command + "\n")
                        # sys.exit(1)
                        no_struct_variants = False #we did not manage to polish the contig, let's try to reassemble it

                    #now retrieve the repolished sequence
                    with open("tmp_repolished.fa", 'r') as repolished :
                        repolished.readline()
                        seq = repolished.readline().strip()


                if not no_struct_variants and s!= 0 and s!= len(names)-1: #let's try to reassemble the reads using neighboring contigs to anchor them
                    
                    print("reassembling ", subcontig, " with ", len(reads[s]), " reads")

                    #now align the reads on the left and right chunks and take the portion of the reads between the two chunks
                    command = "minimap2 -cx map-pb --secondary=no tmp_left.fa tmp_reads.fa > tmp_left.paf 2> trash.txt"
                    minimap = os.system(command)
                    if minimap != 0 :
                        print("Error while running minimap2: " + command + "\n")
                        sys.exit(1)
                    
                    command = "minimap2 -cx map-pb --secondary=no tmp_right.fa tmp_reads.fa > tmp_right.paf 2> trash.txt"
                    minimap = os.system(command)
                    if minimap != 0 :
                        print("Error while running minimap2: " + command + "\n")
                        sys.exit(1)

                    #retrieve the coordinates of the reads mapping on the left chunk
                    left_coordinates = {}
                    with open("tmp_left.paf", 'r') as paf :
                        for line in paf :
                            ls = line.strip().split('\t')
                            if int(ls[11]) == 60 and int(ls[8]) >= int(ls[6])-10: #ls[6] == ls[8] means the read maps to the very end of the contig
                                left_coordinates[ls[0]] = (int(ls[2]), int(ls[3]))

                    #retrieve the coordinates of the reads mapping on the right chunk
                    right_coordinates = {}
                    with open("tmp_right.paf", 'r') as paf :
                        for line in paf :
                            ls = line.strip().split('\t')
                            #if quality of the mapping is good
                            if int(ls[11]) == 60 and int(ls[7]) < 10: #ls[7] == 0 means the read maps to the very beginning of the contig to the right: that's what we want
                                right_coordinates[ls[0]] = (int(ls[2]), int(ls[3]))

                    #now retrieve the reads that are between the two chunks
                    reads_between = {}
                    best_read = "" #that's to measure the read that is best anchored on the sides
                    length_left_and_right = 0
                    idx = 0
                    for read in reads[s] :
                        if read in left_coordinates and read in right_coordinates :
                            # print("read ", read, " is between ", left_coordinates[read], " and ", right_coordinates[read])
                            if left_coordinates[read][0] < right_coordinates[read][0] :
                                reads_between[read] = (max(int(left_coordinates[read][0]), int(left_coordinates[read][1])), 
                                                       min(int(right_coordinates[read][0]), int(right_coordinates[read][1])))
                            else :
                                reads_between[read] = (max(int(right_coordinates[read][0]), int(right_coordinates[read][1])), 
                                                       min(int(left_coordinates[read][0]), int(left_coordinates[read][1])))
                                
                            if reads_between[read][1] - reads_between[read][0] > length_left_and_right :
                                length_left_and_right = reads_between[read][1] - reads_between[read][0]
                                best_read = read
                            idx += 1
                    
                    # print("reads between: ", reads_between)
                    # print("read between: ", [i[1]-i[0] for i in reads_between.values()])

                    #create the list of reads to use for polishing by extracting the reads from the fastq file, cutting them using reads_between and write them to a temporary file
                    f = open("tmp_reads_cut.fa", 'w')
                    
                    f_toPolish = open("tmp_toPolish.fa", 'w')
                    with open(fastq_file, 'r') as fastq :
                        for read in reads_between :
                            fastq.seek(reads_position[read])
                            line = fastq.readline()
                            if read == best_read :
                                contig_seq = line[max(0,reads_between[read][0]-500):min(reads_between[read][1]+500, len(line))]
                                # if orientations_of_reads[read] == "0":
                                #     contig_seq = reverse_complement(contig_seq)
                                f_toPolish.write(">" + read + "\n")
                                f_toPolish.write(contig_seq + "\n") #take a little margin to anchor the contig on both sides
                                
                            else :
                                f.write(">" + read + "\n")
                                f.write(line[max(0,reads_between[read][0]-500):min(reads_between[read][1]+500, len(line))] + "\n")

                    f.close()
                    f_toPolish.close()
                        
                    #now polish f_toPolish with tmp_reads_cut.fa using racon
                    command = "minimap2 -x map-pb -t " + str(threads) + " tmp_toPolish.fa tmp_reads_cut.fa > tmp_toPolish.paf 2> trash.txt"
                    minimap = os.system(command)
                    if minimap != 0 :
                        print("Error while running minimap2: " + command + "\n")
                        sys.exit(1)
                    
                    #check if the alignment is empty
                    empty = True
                    with open("tmp_toPolish.paf", 'r') as paf :
                        for line in paf :
                            ls = line.strip().split('\t')
                            #check if it aligns on more or less the whole read
                            if int(ls[8])-int(ls[7]) > 0.8*int(ls[6]) :
                                empty = False
                                break
                

                    if not empty :
                        command = "racon -w 50 -t " + str(threads) + " tmp_reads_cut.fa tmp_toPolish.paf tmp_toPolish.fa > tmp_repolished.fa 2>trash.txt"
                        racon = os.system(command)
                        if racon != 0 :
                            #polishign failed, fall back sequence
                            seq = contig_seq
                            # print("Error while running racon: " + command + "\n")
                            # sys.exit(1)

                        else:

                            #now retrieve the repolished sequence, realign it one last time against left and right and store it in the segment
                            command = "minimap2 -cx map-pb --secondary=no tmp_left.fa tmp_repolished.fa > tmp_left.paf 2> trash.txt"
                            minimap = os.system(command)
                            if minimap != 0 :
                                print("Error while running minimap2: " + command + "\n")
                                sys.exit(1)
                            command = "minimap2 -cx map-pb --secondary=no tmp_right.fa tmp_repolished.fa > tmp_right.paf 2> trash.txt"
                            minimap = os.system(command)
                            if minimap != 0 :
                                print("Error while running minimap2: " + command + "\n")
                                sys.exit(1)

                            #retrieve the coordinates of the reads mapping on the left chunk
                            reversed_seq = False
                            left_coordinates = (0,0)
                            with open("tmp_left.paf", 'r') as paf :
                                for line in paf :
                                    line = line.strip().split('\t')
                                    left_coordinates = (int(line[2]), int(line[3]))
                                    if line[5] == "-":
                                        reversed_seq = True
                                    break
                            right_coordinates = (0,0)
                            with open("tmp_right.paf", 'r') as paf :
                                for line in paf :
                                    line = line.strip().split('\t')
                                    right_coordinates = (int(line[2]), int(line[3]))
                                    break

                            if reversed_seq :
                                left_coordinates, right_coordinates = right_coordinates, left_coordinates

                            if right_coordinates == (0,0) or left_coordinates == (0,0) :
                                #problem in the polishing
                                print("DEBUG code 3309")
                                seq = None
                            else:
                                #now retrieve the repolished sequence between left_coordinates and right_coordinates
                                with open("tmp_repolished.fa", 'r') as repolished :
                                    repolished.readline()
                                    seq = repolished.readline().strip()
                                    # if left_coordinates == (0,0):
                                    #     left_coordinates = (0, min(500, len(seq)))
                                    # if right_coordinates == (0,0): #should not happen, but could if bad polishing
                                    #     right_coordinates=(max(min(500, len(seq)), len(seq)-500), len(seq))

                                    seq = seq[max(left_coordinates[0], left_coordinates[1])-1:min(right_coordinates[0], right_coordinates[1])+1]
                                    if reversed_seq :
                                        # print("REVERSIIIIIING !", orientations[s])
                                        seq = reverse_complement(seq)
                                        # sys.exit(1)
                                
                                # print("repolished sequence: ", seq)
                    #because we made sure the orientation was positive when choosing left and right
                    segment.set_orientation(s, 1)
                
                if seq is not None :
                    seqs[s] = seq

       

        
        segment.set_sequences(seqs)

    #remove temporary files
    # os.system("rm tmp*")
    return segments






