import re
import sys
import segment as sg
from input_output import read_GAF

class Path :

    def __init__(self, contigs, orientations) :
        if len(contigs) != len(orientations) :
            raise ValueError("ERROR in simple_unzip.py: iiox")
        self.__contigs = contigs
        self.__orientations = []
        for i in orientations :
            if i == ">":
                self.__orientations.append(1)
            else :
                self.__orientations.append(0)

    def __len__(self):
        return len(self.__contigs)
    
    def __str__(self):
        string = ""
        for c in range(len(self.__contigs)) :
            string += "<>"[self.__orientations[c]] + str(self.__contigs[c].names)+ ", "
        return string
    
    def get_contigs(self):
        return self.__contigs
    
    def get_orientations(self):
        return self.__orientations
    
    def replace(self, contig_before, contig_after) :
        for c in range(len(self.__contigs)):
            if self.__contigs[c] == contig_before :
                self.__contigs[c] = contig_after
    
    def cancel(self, contig): #empty the path if it contains this contig
        co = 0
        for c in self.__contigs :
            if c.ID == contig.ID :
                # print("Cancelling ffry ", self, " ", contig.names)
                # if co > 0 and co < len(self.__contigs)-1 :
                #     print("ERROR: cancelddling a path in the middle of it: ", self)
                #     sys.exit()
                self.__contigs = []
                self.__orientations = []
                return
            co += 1

    #check if the contigs are actually linked
    def split_if_invalid(self) :
        all_coherent_subpaths = []
        last_index = 0
        c = 0
        while c < len(self.__contigs)-1:
            index = sg.find_this_link(self.__contigs[c+1], 1-self.__orientations[c+1] , self.__contigs[c].links[self.__orientations[c]], self.__contigs[c].otherEndOfLinks[self.__orientations[c]])
            if index == -1 :
                # print("No link between ", self.__contigs[c+1].names , " and ", self.__contigs[c].names)
                # self.__contigs = []
                # self.__orientations = []
                all_coherent_subpaths += [Path(self.__contigs[last_index:c+1], ["<>"[i] for i in self.__orientations[last_index:c+1]])]
                last_index = c+1
            c+= 1

        all_coherent_subpaths += [Path(self.__contigs[last_index:c+1], ["<>"[i] for i in self.__orientations[last_index:c+1]])]
        return all_coherent_subpaths


def simple_unzip(segments, names, gafFile) :

    lines = []
    print("Reading the gaf file...")
    read_GAF(gafFile, 0, 0, lines)
    print("Finished going through the gaf file.")

    on_which_paths_is_this_contig = {}
    for s in segments :
        on_which_paths_is_this_contig[s] = []

    # translate the paths in the GAF as paths in the graph
    paths = []
    for line in lines :
        #split the line on '>' and '<' to get the path
        cont = re.split('[><]' , line)
        orientations = "".join(re.findall("[<>]", line))
        del cont[0] #because the first element is always ''
        contigs = [segments[names[i]] for i in cont] 
        
        paths.append(Path(contigs, orientations))

        #     sys.exit()

    #make sure all paths are coherent
    new_paths = []
    for p in paths:
        # if "edge_111@0@0_42000_717" in str(p) :
        #     print("Here is sd ddisi Path: ", p)
        new_paths += [i for i in p.split_if_invalid()]
    paths = new_paths

    pa = 0
    for p in paths :
        co = 0
        for c in p.get_contigs() :
            on_which_paths_is_this_contig[c].append((pa, co))
            # if c.names == ['edge_111@0@0_42000_717'] :
            #     print("my litt segment is here! ", p)
            co += 1
        pa += 1


    toDelete = set()
    go_on = True
    while go_on :
        go_on = False
        for segment in segments :

            segment_to_duplicate = False
            #see if it should be duplicated
            if len(segment.links[0]) > 1 or len(segment.links[1]) > 1 : 

                #list the paths going through the contig
                pairs = {}
                pair_to_paths = {}
                for p in on_which_paths_is_this_contig[segment] :
                    # if segment.names == ['edge_111@0@0_42000_717'] :
                    #     print("ciciuddou ", paths[p[0]], " ", p[1], " ", len(path)-1)
                    path = paths[p[0]]
                    if p[1] > 0 and p[1] < len(path)-1 : #we want paths that go through the contigs
                        #find the two links that are supported by this path
                        contigs = path.get_contigs()
                        orientations = path.get_orientations()

                        contig_left = contigs[p[1]-1]
                        end_left = orientations[p[1]-1]
                        index_left = sg.find_this_link(contig_left, end_left, segment.links[1-orientations[p[1]]], segment.otherEndOfLinks[1-orientations[p[1]]])
                        if index_left == -1 :
                            print("ERROR: From semgent ", segment.names, ", did not find ", contig_left.names, " among ", \
                                  [i.names for i in segment.links[1-orientations[p[1]]]], \
                                  " ", path)
                            sys.exit()

                        contig_right = contigs[p[1]+1]
                        end_right = 1-orientations[p[1]+1]
                        index_right = sg.find_this_link(contig_right, end_right, segment.links[orientations[p[1]]], segment.otherEndOfLinks[orientations[p[1]]])

                        pair = (index_left, index_right)
                        if orientations[p[1]] == 0 :
                            pair = (index_right, index_left)
                        if pair not in pairs :
                            pairs[pair] = 0
                            pair_to_paths[pair] = []
                        pairs[pair] += 1
                        pair_to_paths[pair].append(p)

                        # if segment.names == ['edge_111@0@0_42000_717'] :
                        #     print("On contig ", segment.names, ", path ", path, " supports the links towards ", segment.links[1-orientations[p[1]]][index_left].names, \
                        #             " and ", segment.links[orientations[p[1]]][index_right].names)

                # if segment.names == ['edge_111@0@0_82000_1350'] :
                #     print("Looking at segment ", segment.names, " ", pairs, " ", [(segment.links[0][i[0]].names, segment.links[1][i[1]].names) for i in pairs])
                #     # print("qldjl present on ", on_which_paths_is_this_contig[segment])

                #test if the pairs are enough to support a duplication
                new_pairs = {}
                for p in pairs.keys() :
                    if pairs[p] >= 3 :
                        new_pairs[p] = pairs[p]
                pairs = new_pairs

                all_links_left = [False for i in range(len(segment.links[0]))]
                all_links_right = [False for i in range(len(segment.links[1]))]

                for pair in pairs.keys() :
                    all_links_left[pair[0]] = True
                    all_links_right[pair[1]] = True


                segment_to_duplicate = all([i for i in all_links_right+all_links_left])

                if segment_to_duplicate :# and segment.names == ['edge_79@0@0'] :
                    for pair in pairs.keys() :
                        # print("On contig ", segment.names, " a pair is ", segment.links[0][pair[0]].names, " ", segment.links[1][pair[1]].names )

                        #create a new segment
                        new_segment = sg.Segment(segment.names, segment.orientations, segment.lengths, segInsideCIGARs=segment.insideCIGARs)
                        sg.add_link(segment.links[0][pair[0]], segment.otherEndOfLinks[0][pair[0]], new_segment, 0, CIGAR = segment.CIGARs[0][pair[0]])
                        sg.add_link(segment.links[1][pair[1]], segment.otherEndOfLinks[1][pair[1]], new_segment, 1, CIGAR = segment.CIGARs[1][pair[1]])

                        for pa in pair_to_paths[pair] :
                            # if pa[0] == 31624:
                            #     print("1055 is heer e ", paths[pa[0]], " ", segment.ID, " ", new_segment.ID, " ", [i.ID for i in paths[pa[0]].get_contigs()])
                            paths[pa[0]].replace(segment, new_segment)
                        on_which_paths_is_this_contig[new_segment] = pair_to_paths[pair]

                        segments.append(new_segment)

                    pa = 0
                    for path in paths :
                        # if segment.names == ['edge_111@0@0_2000_1358'] :
                        # if pa == 31624 :
                            # print("On contig ", segment.names, " the pairs are ", [(segment.links[0][pair[0]].names, segment.links[1][pair[1]].names) for pair in pairs.keys()] )
                        #     print("path ", path, " ", segment.names, " ", pairs, " ", pa)
                        path.cancel(segment)
                        pa += 1

                    segment.cut_all_links()
                    toDelete.add(segment)
                    on_which_paths_is_this_contig[segment] = []

                    go_on = True
            
    delIdx = 0
    while delIdx < len(segments) :
        if segments[delIdx] in toDelete :
            del segments[delIdx]
            del on_which_paths_is_this_contig[segments[delIdx]]
        else :
            delIdx += 1

    return segments