#!/usr/bin/env python3
# ==============================================================
# Lars Gabriel
#
# Data structure for a graph that detects overlapping transcripts of multiple genome annotations
# ==============================================================

class Node:

    def __init__(self, a_id, t_id):
        self.id = '{}_{}'.format(a_id, t_id)
        self.transcript_id = t_id
        self.anno_id = a_id
        self.edge_to = []

class Graph:

    def __init__(self, genome_anno_lst):
        self.nodes = {}
        self.anno = {}
        self.duplicates = {}
        self.component_list = []
        self.init_anno(genome_anno_lst)

    def init_anno(self, genome_anno_lst):
        #make sure that the genome_anno ids are unique
        counter = 0
        for ga in genome_anno_lst:
            if ga.id in self.anno.keys():
                counter += 1
                new_id = "duplicate_anno_{}".format(counter)
                self.duplicates.update({new_id : ga.id})
                ga.change_id(new_id)
            self.anno.update({ga.id : ga})

    def build(self):
        # build graph
        # put all transcripts in one list and sort list by start coordinates
        transcripts = []
        for k in self.anno.keys():
            transcripts += self.anno[k].get_transcript_list()
        transcripts = sorted(transcripts, key=lambda t:t.start)

        for i in range(0, len(transcripts)):
            # add a node for each transcript t
            key = '{}_{}'.format(transcripts[i].source_anno, transcripts[i].id)
            self.nodes.update({key : Node(transcripts[i].source_anno, \
                transcripts[i].id)})
            # detect overlapping transcripts and add an edge to them
            # find overlapping transcripts t_j with t_j.end <= t.start
            j = i
            while True:
                j -= 1
                if j < 0:
                    break
                if transcripts[j].end < transcripts[i].start:
                    break
                self.nodes[key].edge_to.append('{}_{}'.\
                    format(transcripts[j].source_anno, transcripts[j].id))
            # find overlapping transcripts t_j with t_j.start <= t.end
            j = i
            while True:
                j += 1
                if j == len(transcripts):
                    break
                if transcripts[j].start > transcripts[i].end:
                    break
                self.nodes[key].edge_to.append('{}_{}'.\
                    format(transcripts[j].source_anno, transcripts[j].id))

        for k in self.nodes.keys():
            print(self.nodes[k].id)
            print(self.nodes[k].transcript_id)
            print(self.nodes[k].anno_id)
            print(self.nodes[k].edge_to)
            print('\n')

    def connected_components(self):
        #returns list of connected components
        visited = []
        self.component_list
        for key in list(self.nodes.keys()):
            print(key)
            component = [key]
            if key in visited:
                continue
            visited.append(key)
            not_visited = self.nodes[key].edge_to.copy()
            component += not_visited
            while not_visited:
                next_node = not_visited.pop()
                visited.append(next_node)
                new_nodes = [n for n in self.nodes[next_node].edge_to if n not in component]
                not_visited += new_nodes
                component += new_nodes
            self.component_list.append(component)

        print(self.component_list)
        return self.component_list
