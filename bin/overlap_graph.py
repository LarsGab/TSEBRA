#!/usr/bin/env python3
# ==============================================================
# Lars Gabriel
#
# Graph for transcripts of multiple genome annotations
# It can detect overlapping transcripts
# Add a feature vector to each node
# Compare nodes with the 'decision rule'
# ==============================================================
from features import Edge_features

class Edge:
    def __init__(self, tx1_id, tx2_id):
        self.node1 = tx1_id
        self.node2 = tx2_id
        self.features = {}
        self.node_to_remove = None

class Node:
    def __init__(self, a_id, t_id):
        self.id = '{};{}'.format(a_id, t_id)
        self.transcript_id = t_id
        self.anno_id = a_id
        self.component_id = None

        # dict of edge_ids of edges that are incident
        # self.edge_to[id of incident Node] = edge_id
        self.edge_to = {}

class Graph:
    def __init__(self, genome_anno_lst, anno_pref='braker2'):
        # self.nodes['anno;txid'] = Node(anno, txid)
        self.nodes = {}

        # self.edges['ei'] = Edge()
        self.edges = {}

        # self.anno[annoid] = Anno()
        self.anno = {}

        # preferred annotation source
        self.anno_pref = anno_pref

        # list of connected graph components
        self.component_list = []

        # subset of all transcripts that weren't excluded by the decision rule
        self.decided_graph = []

        # dict of duplicate genome annotation ids to new ids
        self.duplicates = {}

        # init annotations, check fpr duplicate ids
        self.init_anno(genome_anno_lst)

    def init_anno(self, genome_anno_lst):
        # make sure that the genome_anno ids are unique
        counter = 0
        for ga in genome_anno_lst:
            if ga.id in self.anno.keys():
                counter += 1
                new_id = "duplicate.anno.{}".format(counter)
                self.duplicates.update({new_id : ga.id})
                ga.change_id(new_id)
            self.anno.update({ga.id : ga})

    def __tx_from_key__(self, key):
        # get transcript of a node
        anno_id, tx_id = key.split(';')
        return self.anno[anno_id].transcripts[tx_id]

    def build(self):
        # create vertex in graph for each transcript
        # tx_start_end[chr] = [tx_id, coord, id for start or end]
        # for every tx one element for start and one for end
        tx_start_end = {}

        for k in self.anno.keys():
            for tx in self.anno[k].get_transcript_list():
                if tx.chr not in tx_start_end.keys():
                    tx_start_end.update({tx.chr : []})
                key = '{};{}'.format(tx.source_anno, \
                    tx.id)
                self.nodes.update({key : Node(tx.source_anno, \
                    tx.id)})
                tx_start_end[tx.chr].append([key, tx.start, 1])
                tx_start_end[tx.chr].append([key, tx.end, 0])

        # detect overlapping transcripts
        edge_count = 0
        for chr in tx_start_end.keys():
            tx_start_end[chr] = sorted(tx_start_end[chr], key=lambda t:(t[1], t[2]))
            open_intervals = []
            for interval in tx_start_end[chr]:
                if interval[2] == 1:
                    open_intervals.append(interval[0])
                else:
                    open_intervals.remove(interval[0])
                    for match in open_intervals:
                        tx1 = self.__tx_from_key__(interval[0])
                        tx2 = self.__tx_from_key__(match)
                        if self.compare_tx_cds(tx1, tx2):
                            new_edge_key = 'e{}'.format(edge_count)
                            edge_count += 1
                            self.edges.update({new_edge_key : Edge(interval[0], match)})
                            self.nodes[interval[0]].edge_to.update({match : new_edge_key})
                            self.nodes[match].edge_to.update({interval[0] : new_edge_key})

    def compare_tx_cds(self, tx1, tx2):
        # check for overlapping cds in two txs
        coords = []
        coords += tx1.get_cds_coords()
        coords += tx2.get_cds_coords()
        coords = sorted(coords, key=lambda c:c[0])
        for i in range(1, len(coords)):
            if coords[i][0] <= coords[i-1][1]:
                return True
        return False

    def print_nodes(self):
        for k in self.nodes.keys():
            print(self.nodes[k].id)
            print(self.nodes[k].transcript_id)
            print(self.nodes[k].anno_id)
            print(self.nodes[k].edge_to.keys())
            print('\n')

    def connected_components(self):
        # returns list of connected components of the graph and adds component id to vertices
        visited = []
        self.component_list = []
        component_index = 0
        for key in list(self.nodes.keys()):
            component = [key]
            if key in visited:
                continue
            visited.append(key)
            not_visited = list(self.nodes[key].edge_to.keys())
            component += not_visited
            while not_visited:
                next_node = not_visited.pop()
                visited.append(next_node)
                new_nodes = [n for n in self.nodes[next_node].edge_to.keys() if n not in component]
                not_visited += new_nodes
                component += new_nodes
            self.component_list.append(component)
            component_index += 1
            for node in component:
                self.nodes[node].component_id = 'g_{}'.format(component_index)
        return self.component_list

    def add_edge_features(self, hintfiles):
        # hintfiles: list of Hintfiles()
        for key in self.edges.keys():
            tx1 = self.__tx_from_key__(self.edges[key].node1)
            tx2 = self.__tx_from_key__(self.edges[key].node2)

            if tx1.start < tx2.start:
                start = tx1.start
            else:
                start = tx2.start
            if tx1.end > tx2.end:
                end = tx1.end
            else:
                end = tx2.end

            evi = []
            for h_file in hintfiles:
                evi += h_file.hints_in_range(start, end, tx1.chr)

            type_dict = {'intron' : ['intron'], 'start_stop' : \
                ['start_codon', 'stop_codon']}
            for type_key in type_dict.keys():
                tx1_lst = []
                tx2_lst = []
                evi_lst = []
                for type in type_dict[type_key]:
                    if type in tx1.transcript_lines.keys():
                        tx1_lst += tx1.transcript_lines[type]
                    if type in tx2.transcript_lines.keys():
                        tx2_lst += tx2.transcript_lines[type]
                    evi_lst += [e for e in evi if e[2] == type]
                self.edges[key].features.update({type_key : Edge_features(tx1_lst, \
                    tx2_lst, evi_lst).feature_vector})

    def decide_edge(self, edge):
        # decision rule
        # compare the features of the nodes connected by edge
        for i in [0,1]:
            for type in ['intron', 'start_stop']:
                if edge.features[type][i] > edge.features[type][i+2]:
                    return edge.node2
                if edge.features[type][i] < edge.features[type][i+2]:
                    return edge.node1
        if edge.node1.split(';')[0] == self.anno_pref:
            return edge.node2
        elif edge.node2.split(';')[0] == self.anno_pref:
            return edge.node1
        return None

    def decide_component(self, component):
        result = component.copy()
        for node_id in component:
            for e_id in self.nodes[node_id].edge_to.values():
                node_to_remove = self.edges[e_id].node_to_remove
                if node_to_remove:
                    if node_to_remove in result:
                        result.remove(node_to_remove)
        return result

    def decide_graph(self):
        for key in self.edges.keys():
            self.edges[key].node_to_remove = self.decide_edge(self.edges[key])
        self.decided_graph = []
        if not self.component_list:
            self.connected_components()
        for component in self.component_list:
            if len(component) > 1:
                self.decided_graph += self.decide_component(component)
            else:
                self.decided_graph += component

    def get_decided_graph(self):
        # returns the result of decide graph as a dict with
        # result[anno_id] = [[tx_ids, new_gene_id]]
        # a connected component in decided graph is a gene in the output
        if not self.decided_graph:
            self.decide_graph()
        result = {}
        for key in self.anno.keys():
            result.update({key : []})
        for node in self.decided_graph:
            anno_id, tx_id = node.split(';')
            result[anno_id].append([tx_id, self.nodes[node].component_id])
        return result
