#author: Roland Faure
#program which bluntifies a GFA 

version = "0.1"
date = "2024-08-22"

import argparse
import sys
import re

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement[base] for base in reversed(seq))

def basic_overlap_removal(gfa_in, gfa_out, trust):
    #parse the graph
    list_of_contigs= {} #each contig is a name: (linksleft, linksright)
    length_of_contigs = {} #associates to each contig its length
    location_of_contigs_in_gfa = {} #associates to each contig its location in the gfa file, in order to quickly access the contigs when recomputing the overlaps
    list_of_links = [] #each link is a tuple (name1, end1, overlap1, name2, end2, overlap2)
    set_of_already_appended_links = set() #to avoid duplicates

    with open(gfa_in, "r") as gfa_open:
        line = gfa_open.readline()
        while line:
            if line.startswith("S"):
                fields = line.strip().split("\t")
                name = fields[1]
                list_of_contigs[name] = ([],[])
                length = len(fields[2])
                location_of_contigs_in_gfa[name] = gfa_open.tell() - len(line)
                length_of_contigs[name] = length
            elif line.startswith("L"):
                fields = line.strip().split("\t")
                name1 = fields[1]
                end1 = int(fields[2]=="+")
                name2 = fields[3]
                end2 = int(fields[4]=="-")
                if (name1, end1, name2, end2) in set_of_already_appended_links:
                    line = gfa_open.readline()
                    continue
                else:
                    set_of_already_appended_links.add((name1, end1, name2, end2))
                    set_of_already_appended_links.add((name2, end2, name1, end1))

                #compute the overlap between the contigs
                cigar = fields[5]

                #compute the length of the overlaps on both contigs
                overlap1 = 0
                overlap2 = 0
                len_string = ""
                for c in cigar:
                    if c.isdigit():
                        len_string += c
                    else:
                        if c == "M":
                            overlap1 += int(len_string)
                            overlap2 += int(len_string)
                        elif c == "I":
                            print("ERROR: bluntify only works with perfect overlaps")
                            overlap2 += int(len_string)
                        elif c == "D":
                            overlap1 += int(len_string)
                        len_string = ""

                if not trust: #then we need to recompute the overlaps
                    #align the two contigs - first fetch the sequences
                    open2 = open(gfa_file, "r")
                    open2.seek(location_of_contigs_in_gfa[name1])
                    seq1 = open2.readline().strip().split("\t")[2]
                    if end1 == 0:
                        seq1_cut = reverse_complement(seq1[:min(int(1.1*overlap1)+2, len(seq1))])#take a bit of margin to be sure to have the full overlap
                    else :
                        seq1_cut = seq1[max(int(-1.1*overlap1)-2, -len(seq1)):]
                    open2.seek(location_of_contigs_in_gfa[name2])
                    seq2 = open2.readline().strip().split("\t")[2]
                    if end2 == 1:
                        seq2_cut = reverse_complement(seq2)[:min(len(seq2), int(1.1*overlap2)+2)]
                    else:
                        seq2_cut = seq2[:min(len(seq2), int(1.1*overlap2)+2)]
                    open2.close()

                    #find the end of the overlaps
                    alignment = edlib.align(seq1_cut, seq2_cut, task="path", mode="SHW")
                    new_overlap2 = alignment["locations"][0][1] + 1
                    alignment= edlib.align(reverse_complement(seq2_cut), reverse_complement(seq1_cut), task="path", mode="SHW")
                    new_overlap1 = alignment["locations"][0][1] + 1

                    print("\n",(name1, end1, overlap1, name2, end2, overlap2))
                    #update the overlaps
                    if overlap1 != None and overlap2 != None:
                        overlap1 = new_overlap1
                        overlap2 = new_overlap2


                    print((name1, end1, overlap1, name2, end2, overlap2))

                list_of_links.append([name1, end1, overlap1, name2, end2, overlap2])
                list_of_contigs[name1][end1].append(len(list_of_links)-1)
                list_of_contigs[name2][end2].append(len(list_of_links)-1)

            line = gfa_open.readline()

    #greedily delete overlaps
    trimmed_lengths = {}#associate to each contig a trimmed length left and right
    for contig in list_of_contigs.keys():
        #look at all the links left and right
        min_length_left = length_of_contigs[contig]
        max_length_left = 0
        if len(list_of_contigs[contig][0])==0:
            min_length_left = 0
        min_length_right = length_of_contigs[contig]
        max_length_right= 0
        if len(list_of_contigs[contig][1]) == 0:
            min_length_right = 0

        for link_idx in list_of_contigs[contig][0]:
            overlap_length = list_of_links[link_idx][2]
            if overlap_length > max_length_left :
                max_length_left = overlap_length
            if overlap_length < min_length_left :
                min_length_left = overlap_length

        for link_idx in list_of_contigs[contig][1]:
            overlap_length = list_of_links[link_idx][2]
            if overlap_length > max_length_right :
                max_length_right = overlap_length
            if overlap_length < min_length_right :
                min_length_right = overlap_length

        trim_left = min(min_length_left, length_of_contigs[contig] - max_length_right)
        trim_right = min(min_length_right ,  length_of_contigs[contig] - max_length_left)

        #update all the links
        for link_idx in list_of_contigs[contig][0]:
            list_of_links[link_idx][2] -= trim_left
            list_of_links[link_idx][5] -= trim_left
        for link_idx in list_of_contigs[contig][1]:
            if list_of_links[link_idx][0] != list_of_links[link_idx][3]:#for self_loops
                list_of_links[link_idx][2] -= trim_right

        trimmed_lengths[contig] = (trim_left, trim_right)

    #now that all contigs are trimmed and the length of overlaps diminished, output the new GFA
    fo = open(gfa_out, "w")
    fi = open(gfa_in, "r")
    for line in fi:
        if line[0] == "S" :
            ls = line.strip().split('\t')
            fo.write("S\t"+ls[1]+"\t"+ls[2][trimmed_lengths[ls[1]][0]:length_of_contigs[ls[1]]-trimmed_lengths[ls[1]][1]]+"\t"+"\t".join(ls[3:])+"\n")
    for contig in list_of_contigs.keys():
        #output all links then delete them 
        for link_idx in list_of_contigs[contig][0]:
            if list_of_links[link_idx][0] != "None" :
                
                fo.write("L\t"+ list_of_links[link_idx][0]+"\t"+ "-+"[list_of_links[link_idx][1]] + "\t" + list_of_links[link_idx][3] + "\t"+\
                         "+-"[list_of_links[link_idx][4]]+ "\t"+str(list_of_links[link_idx][2])+"M\n")
                list_of_links[link_idx][0] = "None"

        for link_idx in list_of_contigs[contig][1]:
            if list_of_links[link_idx][0] != "None" :
                
                fo.write("L\t"+ list_of_links[link_idx][0]+"\t"+ "-+"[list_of_links[link_idx][1]] + "\t" + list_of_links[link_idx][3] + "\t"+\
                         "+-"[list_of_links[link_idx][4]]+ "\t"+str(list_of_links[link_idx][2])+"M\n")
                list_of_links[link_idx][0] = "None"

    fo.close()
    fi.close()

def fancier_overlap_removal(gfa_in, gfa_out, say_0M_no_matter_what = False, short_contig_length = 0):
    #parse the graph
    list_of_contigs= {} #each contig is a name: (linksleft, linksright)
    length_of_contigs = {} #associates to each contig its length
    list_of_links = [] #each link is a tuple (name1, end1, overlap1, name2, end2, overlap2)
    set_of_already_appended_links = set() #to avoid duplicates
    location_of_contigs_in_gfa = {}
    coverage = {} #associates to each contig the coverage

    with open(gfa_in, "r") as gfa_open:
        line = gfa_open.readline()
        while line:
            if line.startswith("S"):
                fields = line.strip().split("\t")
                name = fields[1]
                list_of_contigs[name] = ([],[])
                length = len(fields[2])
                length_of_contigs[name] = length
                location_of_contigs_in_gfa[name] = gfa_open.tell() - len(line)
                #find the coverage field
                for field in fields:
                    if field.upper().startswith("DP:"):
                        coverage[name] = field
                        break
            elif line.startswith("L"):
                fields = line.strip().split("\t")
                name1 = fields[1]
                end1 = int(fields[2]=="+")
                name2 = fields[3]
                end2 = int(fields[4]=="-")
                if (name1, end1, name2, end2) in set_of_already_appended_links:
                    line = gfa_open.readline()
                    continue
                else:
                    set_of_already_appended_links.add((name1, end1, name2, end2))
                    set_of_already_appended_links.add((name2, end2, name1, end1))

                #compute the overlap between the contigs
                cigar = fields[5]

                #compute the length of the overlaps on both contigs
                overlap1 = 0
                overlap2 = 0
                len_string = ""
                for c in cigar:
                    if c.isdigit():
                        len_string += c
                    else:
                        if c == "M":
                            overlap1 += int(len_string)
                            overlap2 += int(len_string)
                        elif c == "I":
                            overlap2 += int(len_string)
                        elif c == "D":
                            overlap1 += int(len_string)
                        len_string = ""

                list_of_links.append([name1, end1, overlap1, name2, end2, overlap2])
                list_of_contigs[name1][end1].append(len(list_of_links)-1)
                list_of_contigs[name2][end2].append(len(list_of_links)-1)

            line = gfa_open.readline()

    new_list_of_contigs = {}
    new_list_of_links = []
    contigs_already_dealt_with = set()
    old_contigs_to_new_contigs = {} #associates to each old contig the [newcontignametotheleft, newcontignametotheright]
    for contig in list_of_contigs.keys():

        #look at all the links left and right
        max_length_left = 0
        max_length_right= 0
        list_breakpoints_left = [] #list of (breakpoint, link)
        list_breakpoints_right = [] #list of (breakpoint, link)
        for link_idx in list_of_contigs[contig][0]:
            overlap_length = list_of_links[link_idx][2]
            if overlap_length > max_length_left :
                max_length_left = overlap_length
            list_breakpoints_left.append((overlap_length, link_idx))
        if (len(list_breakpoints_left)) == 0:
            list_breakpoints_left = [(0,-1)]

        for link_idx in list_of_contigs[contig][1]:
            overlap_length = list_of_links[link_idx][2]
            if overlap_length > max_length_right :
                max_length_right = overlap_length
            list_breakpoints_right.append((length_of_contigs[contig]-overlap_length, link_idx))
        if (len(list_breakpoints_right))== 0:
            list_breakpoints_right = [(length_of_contigs[contig],-1)]

        if max_length_left+max_length_right < length_of_contigs[contig]: #then there is no risks of creating new contigs
            contigs_already_dealt_with.add(contig)

            breakpoints = list_breakpoints_left+list_breakpoints_right
            breakpoints.sort(key= lambda x : x[0]) #sort by breakpoints, then put the links left berfore the links right
            
            # create all the new contigs
            new_contig_name = ""
            for sub in range(0,len(breakpoints)):  
                    
                if sub > 0 and breakpoints[sub-1][0] != breakpoints[sub][0] :

                    new_contig_name = contig + "_" + str(breakpoints[sub-1][0]) + "_" + str(breakpoints[sub][0])
                    n=2
                    while sub-n >= 0 and breakpoints[sub-1][0] == breakpoints[sub-n][0]:
                        n+=1
                    if sub-n >= 0:
                        past_contig_name = contig + "_" + str(breakpoints[sub-n][0]) + "_" + str(breakpoints[sub-1][0])

                        new_list_of_links.append([new_contig_name, 0, 0, past_contig_name, 1, 0])
                        new_list_of_contigs[new_contig_name] = [[len(new_list_of_links)-1],[]]
                    else:
                        new_list_of_contigs[new_contig_name] = [[],[]]
                
                if sub == 0:
                    n = 1
                    while str(breakpoints[sub][0]) == str(breakpoints[sub+n][0]):
                        n+=1
                    future_contig_name = contig + "_" + str(breakpoints[sub][0]) + "_" + str(breakpoints[sub+n][0])
                    old_contigs_to_new_contigs[contig] = [future_contig_name, ""]
                    # print("new contig name ", future_contig_name)
                if sub == len(breakpoints)-1:
                    old_contigs_to_new_contigs[contig][1] = new_contig_name

            #add the links
            for sub in range(0,len(breakpoints)):
                # if contig == "utg000037l":
                #     print("adding link ", list_of_links[breakpoints[sub][1]])

                #add the links, with placeholders for the names we don't know yet
                if breakpoints[sub][1] != -1:
                    link = list_of_links[breakpoints[sub][1]].copy()
                    link[2] = 0
                    link[5] = 0
                    if link[0] == -1:
                        continue

                    if (link[0] == contig and link[1] == 0) or (link[3] == contig and link[4] == 0):
                        n = 1
                        while str(breakpoints[sub][0]) == str(breakpoints[sub+n][0]):
                            n+=1
                        future_contig_name = contig + "_" + str(breakpoints[sub][0]) + "_" + str(breakpoints[sub+n][0])
                        if link[0] == contig and link[1] == 0:
                            link[0] = future_contig_name
                        elif link[3] == contig and link[4] == 0:
                            link[3] = future_contig_name
                        new_list_of_links.append(link)
                        new_list_of_contigs[future_contig_name][0].append(len(new_list_of_links)-1)
                        
                    else:
                        n=1
                        while breakpoints[sub][0] == breakpoints[sub-n][0]:
                            n+=1
                        new_contig_name = contig + "_" + str(breakpoints[sub-n][0]) + "_" + str(breakpoints[sub][0])
                        if link[0] == contig and link[1] == 1:
                            link[0] = new_contig_name
                        if link[3] == contig and link[4] == 1:
                            link[3] = new_contig_name
                        new_list_of_links.append(link)
                        new_list_of_contigs[new_contig_name][1].append(len(new_list_of_links)-1)
    
                    # if "utg000037l" in link[0] or "utg000037l" in link[3]:
                    #     # print(breakpoints)
                    #     print("Here is a link I want to append, ", contig, " ", link, " original: ", list_of_links[breakpoints[sub][1]])


                    # #delete the overlap in the link, to avoid redetaching it
                    # list_of_links[breakpoints[sub][1]][2] = 0
                    # list_of_links[breakpoints[sub][1]][5] = 0
                    #delete the link, to avoid reusing it
                    list_of_links[breakpoints[sub][1]][0] = -1
                    list_of_links[breakpoints[sub][1]][3] = -1

        else: #overlap everywhere, just copy the contig
            #copy the contig
            new_name = contig + "_0_" + str(length_of_contigs[contig])
            new_list_of_contigs[new_name] = ([],[])
            old_contigs_to_new_contigs[contig] = [new_name, new_name]
            # #copy the links left
            # for link_idx in list_of_contigs[contig][0]:
            #     new_list_of_links.append(list_of_links[link_idx].copy())
            #     new_list_of_links[-1][0] = new_name
            #     new_list_of_contigs[new_name][0].append(len(new_list_of_links)-1)
            # #copy the links right
            # for link_idx in list_of_contigs[contig][1]:
            #     new_list_of_links.append(list_of_links[link_idx].copy())
            #     new_list_of_links[-1][3] = new_name
            #     new_list_of_contigs[new_name][1].append(len(new_list_of_links)-1)

    # print(old_contigs_to_new_contigs)
    # go through all the links and update them with their new names
    for link in new_list_of_links:
        if link[0] in old_contigs_to_new_contigs:
            if link[1] == 0:
                link[0] = old_contigs_to_new_contigs[link[0]][0]
            else:
                link[0] = old_contigs_to_new_contigs[link[0]][1]
        if link[3] in old_contigs_to_new_contigs:
            if link[4] == 0:
                link[3] = old_contigs_to_new_contigs[link[3]][0]
            else:
                link[3] = old_contigs_to_new_contigs[link[3]][1]

    #now output the new GFA
    fo = open(gfa_out, "w")
    new_contigs_length = {}
    for contig in new_list_of_contigs.keys():
        #get the sequence of the full contig
        open2 = open(gfa_in, "r")
        original_name = "_".join(contig.split("_")[:-2])
        open2.seek(location_of_contigs_in_gfa[original_name])
        seq = open2.readline().strip().split("\t")[2]
        open2.close()

        #now take the last two parts of the name to get the breakpoints
        left = contig.split("_")[-2]
        right = contig.split("_")[-1]

        #output the new contig
        if int(right) - int(left) >= short_contig_length:
            fo.write("S\t"+contig+"\t"+seq[int(left):int(right)]+"\t")
            if original_name in coverage:
                fo.write(coverage[original_name]+"\n")
            else:
                fo.write("\n")
        new_contigs_length[contig] = int(right) - int(left)

    for contig in new_list_of_contigs.keys():
        #output all links then delete them 
        for link_idx in new_list_of_contigs[contig][0]:
            if new_list_of_links[link_idx][0] != "None" :
                if not say_0M_no_matter_what or new_list_of_links[link_idx][2] == 0:
                    fo.write("L\t"+ new_list_of_links[link_idx][0]+"\t"+ "-+"[new_list_of_links[link_idx][1]] + "\t" + new_list_of_links[link_idx][3] + "\t"+\
                        "+-"[new_list_of_links[link_idx][4]]+ "\t"+str(new_list_of_links[link_idx][2])+"M\n")
                new_list_of_links[link_idx][0] = "None"

        for link_idx in new_list_of_contigs[contig][1]:
            if new_list_of_links[link_idx][0] != "None" :
                if not say_0M_no_matter_what or new_list_of_links[link_idx][2] == 0:
                    fo.write("L\t"+ new_list_of_links[link_idx][0]+"\t"+ "-+"[new_list_of_links[link_idx][1]] + "\t" + new_list_of_links[link_idx][3] + "\t"+\
                        "+-"[new_list_of_links[link_idx][4]]+ "\t"+str(new_list_of_links[link_idx][2])+"M\n")
                new_list_of_links[link_idx][0] = "None"

    fo.close()
                
def main():
    parser = argparse.ArgumentParser(description="Bluntify a GFA file")
    parser.add_argument("input", help="input GFA file")
    parser.add_argument("output", help="output GFA file")
    parser.add_argument("-t", "--trust", action="store_true", help="trust the overlaps")
    parser.add_argument("-n", "--no_overlaps", action="store_true", help="At the end, cut all links which we have not been able to bluntify")
    parser.add_argument("--version", action="version", version="%(prog)s " + version)
    args = parser.parse_args()

    if not args.trust:
        import edlib #for the alignment


    intermediate_gfa = "intermediate_gfa.tmp.gfa"
    basic_overlap_removal(args.input, intermediate_gfa, args.trust)
    fancier_overlap_removal(intermediate_gfa, args.output, say_0M_no_matter_what=args.no_overlaps)

    
if __name__ == "__main__":
    main()
    
        

